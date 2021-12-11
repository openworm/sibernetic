/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2021 OpenWorm.
 * http://openworm.org
 *
 * Substantial portions of alloc_picture, open_video, add_stream,
 * encoder_start, encoder_finish, get_video_frame, write_frame
 * Copyright (c) 2003 Fabrice Bellard
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

#if FFMPEG
#include <cstddef>
#include <cstdio>
#include <iostream>

#include "owVideoWriter.h"

owVideoWriter::owVideoWriter(
        const std::string & formatName,
        const std::string & fileName,
        int fps,
        int width,
        int height) {
    this->width = width;
    this->height = height;
    this->fps = fps;

    encoder_start(
            fileName.c_str(),
            formatName == ""? NULL:formatName.c_str());
}

AVFrame *owVideoWriter::alloc_picture(enum AVPixelFormat pix_fmt, int width, int height) {
    AVFrame *picture;
    int ret;

    picture = av_frame_alloc();
    if (!picture)
        return NULL;

    picture->format = pix_fmt;
    picture->width  = width;
    picture->height = height;

    /* allocate the buffers for the frame data */
    ret = av_frame_get_buffer(picture, 0);
    if (ret < 0) {
        std::cerr << "Could not allocate frame data." << std::endl;
        exit(1);
    }

    return picture;
}

void owVideoWriter::open_video(AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg) {
    int ret;
    AVCodecContext *c = ost->enc;

    /* open the codec */
    ret = avcodec_open2(c, codec, NULL);

    if (ret < 0) {
        std::cerr << "Could not open video codec: " << ERRMSG_AV(ret) << std::endl;
        exit(1);
    }

    /* allocate and init a re-usable frame */
    ost->frame = alloc_picture(c->pix_fmt, c->width, c->height);
    if (!ost->frame) {
        std::cerr << "Could not allocate video frame" << std::endl;
        exit(1);
    }

    /* copy the stream parameters to the muxer */
    ret = avcodec_parameters_from_context(ost->st->codecpar, c);
    if (ret < 0) {
        std::cerr << "Could not copy the stream parameters" << std::endl;
        exit(1);
    }
}

void owVideoWriter::add_stream(OutputStream *ost,
                       AVCodec **codec,
                       enum AVCodecID codec_id) {
    AVCodecContext *c;

    /* find the encoder */
    *codec = avcodec_find_encoder(codec_id);
    if (!(*codec)) {
        std::cerr << "Could not find encoder for '" << avcodec_get_name(codec_id) << "'" << std::endl;
        exit(1);
    }

    ost->st = avformat_new_stream(oc, NULL);
    if (!ost->st) {
        std::cerr << "Could not allocate stream" << std::endl;
        exit(1);
    }
    ost->st->id = oc->nb_streams - 1;
    c = avcodec_alloc_context3(*codec);
    if (!c) {
        std::cerr << "Could not alloc an encoding context" << std::endl;
        exit(1);
    }
    ost->enc = c;

    switch ((*codec)->type) {
    case AVMEDIA_TYPE_VIDEO:
        c->codec_id = codec_id;

        c->bit_rate = 400000;
        c->width = width;
        c->height = height;
        c->gop_size = 10;
        c->max_b_frames = 1;
        c->pix_fmt = AV_PIX_FMT_YUV420P;

        /* timebase: This is the fundamental unit of time (in seconds) in terms
         * of which frame timestamps are represented. For fixed-fps content,
         * timebase should be 1/framerate and timestamp increments should be
         * identical to 1. */
        ost->st->time_base = (AVRational){ 1, fps };
        c->time_base       = ost->st->time_base;

        if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
            /* Needed to avoid using macroblocks in which some coeffs overflow.
             * This does not happen with normal video, it just happens here as
             * the motion of the chroma plane does not match the luma plane. */
            c->mb_decision = 2;
        }
    break;

    default:
        break;
    }

    /* Some formats want stream headers to be separate. */
    if (oc->oformat->flags & AVFMT_GLOBALHEADER) {
        c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
    }
}

void owVideoWriter::encoder_start(const char *filename, const char *codec_id) {
    AVCodec *codec;
    int ret;
#if ( LIBAVFORMAT_VERSION_INT < AV_VERSION_INT(58,9,100) )
    av_register_all();
#endif

    avformat_alloc_output_context2(&oc, NULL, codec_id, filename);
    if (!oc) {
        std::cerr << "Output format " << codec_id << " not found for " << filename << std::endl;
        const AVOutputFormat *item = NULL;
        void *iter = NULL;
        std::cerr << "Available formats: " << std::endl;
        while ((item = av_muxer_iterate(&iter)) != NULL) {
            std::cerr << item->name << " : " << item->long_name
                << " : " << item->extensions << std::endl;
        }
    }

    fmt = oc->oformat;

    add_stream(&video_st, &codec, fmt->video_codec);
    open_video(codec, &video_st, NULL);
    /* open the output file, if needed */
    if (!(fmt->flags & AVFMT_NOFILE)) {
        ret = avio_open(&oc->pb, filename, AVIO_FLAG_WRITE);
        if (ret < 0) {
            std::cerr << "Could not open '" << filename << "': " << ERRMSG_AV(ret) << std::endl;
            exit(1);
        }
    }

    /* Write the stream header, if any. */
    ret = avformat_write_header(oc, NULL);
    if (ret < 0) {
        std::cerr << "Error occurred when opening output file: " << ERRMSG_AV(ret) << std::endl;
        exit(1);
    }
}

void owVideoWriter::close_stream() {
    OutputStream *ost = &video_st;
    avcodec_free_context(&ost->enc);
    av_frame_free(&ost->frame);
    sws_freeContext(ost->sws_ctx);
}

void owVideoWriter::encoder_finish() {
    /* Write the trailer, if any. The trailer must be written before you
     * close the CodecContexts open when you wrote the header; otherwise
     * av_write_trailer() may try to use memory that was freed on
     * av_codec_close(). */
    av_write_trailer(oc);

    /* Close each codec. */
    close_stream();

    if (!(fmt->flags & AVFMT_NOFILE)) {
        /* Close the output file. */
        avio_closep(&oc->pb);
    }

    /* free the stream */
    avformat_free_context(oc);
}

void owVideoWriter::glread_rgb() {
    size_t cur_gl, cur_rgb, nvals;
    const size_t format_nchannels = 4;
    nvals = format_nchannels * width * height;
    pixels = (GLubyte*) realloc((void*)pixels, nvals * sizeof(GLubyte));
    rgb = (uint8_t*) realloc((void*)rgb, nvals * sizeof(uint8_t));
    /* Get RGBA to align to 32 bits instead of just 24 for RGB. May be faster for FFmpeg. */
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            cur_gl  = format_nchannels * (width * (height - i - 1) + j);
            cur_rgb = format_nchannels * (width * i + j);
            for (size_t k = 0; k < format_nchannels; k++) {
                rgb[cur_rgb + k] = pixels[cur_gl + k];
            }
        }
    }
}

AVFrame *owVideoWriter::get_video_frame() {
    OutputStream *ost = &video_st;
    AVCodecContext *c = ost->enc;
    const int in_linesize[1] = { 4 * c->width };

    /* when we pass a frame to the encoder, it may keep a reference to it
     * internally; make sure we do not overwrite it here */
    int ret = av_frame_make_writable(ost->frame);
    if (ret < 0) {
        std::cerr << "Unable to make the output frame writable: " << ERRMSG_AV(ret) << std::endl;
        exit(1);
    }

    /* as we only generate a YUV420P picture, we must convert it
     * to the codec pixel format if needed */
    if (!ost->sws_ctx) {
        ost->sws_ctx = sws_getCachedContext(ost->sws_ctx,
                c->width, c->height, AV_PIX_FMT_RGB32,
                c->width, c->height, AV_PIX_FMT_YUV420P,
                0, NULL, NULL, NULL);
        if (!ost->sws_ctx) {
            std::cerr << "Could not initialize the conversion context" << std::endl;
            exit(1);
        }
    }
    glread_rgb();
    sws_scale(ost->sws_ctx, (const uint8_t * const *) &rgb,
              in_linesize, 0, c->height, ost->frame->data,
              ost->frame->linesize);

    ost->frame->pts = ost->next_pts++;

    return ost->frame;
}

int owVideoWriter::write_frame(AVCodecContext *c,
                       AVStream *st, AVFrame *frame) {
    int ret;

    // send the frame to the encoder
    ret = avcodec_send_frame(c, frame);
    if (ret < 0) {
        std::cerr << "Error sending a frame to the encoder: " << ERRMSG_AV(ret) << std::endl,
        exit(1);
    }

    while (ret >= 0) {
        AVPacket pkt = { 0 };

        ret = avcodec_receive_packet(c, &pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF) {
            break;
        } else if (ret < 0) {
            std::cerr << "Error encoding a frame: " << ERRMSG_AV(ret) << std::endl;
            exit(1);
        }

        /* rescale output packet timestamp values from codec to stream timebase */
        av_packet_rescale_ts(&pkt, c->time_base, st->time_base);
        pkt.stream_index = st->index;

        /* Write the compressed frame to the media file. */
        ret = av_interleaved_write_frame(oc, &pkt);
        av_packet_unref(&pkt);
        if (ret < 0) {
            std::cerr << "Error while writing output packet: " << ERRMSG_AV(ret) << std::endl;
            exit(1);
        }
    }

    return ret == AVERROR_EOF ? 1 : 0;
}

int owVideoWriter::write_video_frame() {
    return write_frame(video_st.enc, video_st.st, get_video_frame());
}

owVideoWriter::~owVideoWriter() {
    encoder_finish();
    free(pixels);
    free(rgb);
}
#endif
