/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2021 OpenWorm.
 * http://openworm.org
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

#ifndef OWVIDEOWRITER_H
#define OWVIDEOWRITER_H
#if FFMPEG

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libavutil/opt.h>
#include <libswscale/swscale.h>
}

#include "owGL.h"

/* Adapted from: ffmpeg example, muxing.c
 * and from
 * https://github.com/cirosantilli/cpp-cheat/blob/19044698f91fefa9cb75328c44f7a487d336b541/ffmpeg/encode.c
 */
typedef struct OutputStream {
    AVStream *st;
    AVCodecContext *enc;

    /* pts of the next frame that will be generated */
    int64_t next_pts;

    AVFrame *frame;

    struct SwsContext *sws_ctx;
} OutputStream;

#define ERRMSG_AV(__retcode) av_make_error_string(averror_strbuf, sizeof(averror_strbuf), (__retcode))


class owVideoWriter {
public:
    /** Starts up the writer, opens a video stream */
    owVideoWriter(const std::string & formatName, const std::string & fileName, int fps, int width, int height);
    /** Save a video frame from the opengl viewport */
    int write_video_frame();
    ~owVideoWriter();

private:
    GLubyte *pixels = NULL;
    AVFormatContext *oc = NULL;
    AVOutputFormat *fmt = NULL;
    OutputStream video_st = { 0 };

    /* Image width */
    int width;
    /* Image height */
    int height;
    /* Video frames per second */
    int fps;

    uint8_t *rgb = NULL;
    char averror_strbuf[AV_ERROR_MAX_STRING_SIZE];
    AVFrame *alloc_picture(enum AVPixelFormat pix_fmt, int width, int height);
    void open_video(AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg);
    int write_frame(AVCodecContext *c, AVStream *st, AVFrame *frame);
    void encoder_start(const char *filename, const char *codec_id);
    AVFrame *get_video_frame();
    void glread_rgb();
    void add_stream(OutputStream *ost,
                           AVCodec **codec,
                           enum AVCodecID codec_id);
    void encoder_finish();
    void close_stream();
};
#endif
#endif
