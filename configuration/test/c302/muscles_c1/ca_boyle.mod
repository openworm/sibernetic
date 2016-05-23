TITLE Mod file for component: Component(id=ca_boyle type=ionChannelHH)

COMMENT

    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.4.4
         org.neuroml.model   v1.4.4
         jLEMS               v0.9.8.4

ENDCOMMENT

NEURON {
    SUFFIX ca_boyle
    USEION ca READ cai,cao WRITE ica VALENCE 2 ? Assuming valence = 2 (Ca ion); TODO check this!!
    
    RANGE gion                           
    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!
    RANGE conductance                       : parameter
    
    RANGE g                                 : exposure
    
    RANGE fopen                             : exposure
    RANGE e_instances                       : parameter
    
    RANGE e_tau                             : exposure
    
    RANGE e_inf                             : exposure
    
    RANGE e_rateScale                       : exposure
    
    RANGE e_fcond                           : exposure
    RANGE e_timeCourse_tau                  : parameter
    
    RANGE e_timeCourse_t                    : exposure
    RANGE e_steadyState_rate                : parameter
    RANGE e_steadyState_midpoint            : parameter
    RANGE e_steadyState_scale               : parameter
    
    RANGE e_steadyState_x                   : exposure
    RANGE f_instances                       : parameter
    
    RANGE f_tau                             : exposure
    
    RANGE f_inf                             : exposure
    
    RANGE f_rateScale                       : exposure
    
    RANGE f_fcond                           : exposure
    RANGE f_timeCourse_tau                  : parameter
    
    RANGE f_timeCourse_t                    : exposure
    RANGE f_steadyState_rate                : parameter
    RANGE f_steadyState_midpoint            : parameter
    RANGE f_steadyState_scale               : parameter
    
    RANGE f_steadyState_x                   : exposure
    RANGE h_alpha                           : parameter
    RANGE h_k                               : parameter
    RANGE h_ca_half                         : parameter
    RANGE h_instances                       : parameter
    RANGE h_SEC                             : parameter
    
    RANGE h_tau                             : exposure
    
    RANGE h_inf                             : exposure
    
    RANGE h_rateScale                       : exposure
    
    RANGE h_fcond                           : exposure
    
    RANGE h_q                               : exposure
    RANGE e_tauUnscaled                     : derived variable
    RANGE f_tauUnscaled                     : derived variable
    RANGE conductanceScale                  : derived variable
    RANGE fopen0                            : derived variable
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!
    
    conductance = 1.0E-5 (uS)
    e_instances = 2 
    e_timeCourse_tau = 0.100026995 (ms)
    e_steadyState_rate = 1 
    e_steadyState_midpoint = -3.3567998 (mV)
    e_steadyState_scale = 6.74821 (mV)
    f_instances = 1 
    f_timeCourse_tau = 150.87999 (ms)
    f_steadyState_rate = 1 
    f_steadyState_midpoint = 25.1815 (mV)
    f_steadyState_scale = -5.0317597 (mV)
    h_alpha = 0.282473 
    h_k = -1.00056E-8 (mM)
    h_ca_half = 6.41889E-8 (mM)
    h_instances = 1 
    h_SEC = 1000 (ms)
}

ASSIGNED {
    
    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel
    v (mV)
    celsius (degC)
    temperature (K)
    eca (mV)
    ica (mA/cm2)
    
    cai (mM)
    
    cao (mM)
    
    
    e_timeCourse_t (ms)                    : derived variable
    
    e_steadyState_x                        : derived variable
    
    e_rateScale                            : derived variable
    
    e_fcond                                : derived variable
    
    e_inf                                  : derived variable
    
    e_tauUnscaled (ms)                     : derived variable
    
    e_tau (ms)                             : derived variable
    
    f_timeCourse_t (ms)                    : derived variable
    
    f_steadyState_x                        : derived variable
    
    f_rateScale                            : derived variable
    
    f_fcond                                : derived variable
    
    f_inf                                  : derived variable
    
    f_tauUnscaled (ms)                     : derived variable
    
    f_tau (ms)                             : derived variable
    
    h_rateScale                            : derived variable
    
    h_inf                                  : derived variable
    
    h_tau (ms)                             : derived variable
    
    h_q                                    : derived variable
    
    h_fcond                                : derived variable
    
    conductanceScale                       : derived variable
    
    fopen0                                 : derived variable
    
    fopen                                  : derived variable
    
    g (uS)                                 : derived variable
    rate_e_q (/ms)
    rate_f_q (/ms)
    
}

STATE {
    e_q  
    f_q  
    
}

INITIAL {
    eca = 40.0
    
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
    e_q = e_inf
    
    f_q = f_inf
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=ca_boyle type=ionChannelHH), from conductanceScaling; null
    ? Path not present in component, using factor: 1
    
    conductanceScale = 1 
    
    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=ca_boyle type=ionChannelHH), from gates; Component(id=e type=gateHHtauInf)
    ? multiply applied to all instances of fcond in: <gates> ([Component(id=e type=gateHHtauInf), Component(id=f type=gateHHtauInf), Component(id=h type=customHGate)]))
    fopen0 = e_fcond * f_fcond * h_fcond ? path based
    
    fopen = conductanceScale  *  fopen0 ? evaluable
    g = conductance  *  fopen ? evaluable
    gion = gmax * fopen 
    
    ica = gion * (v - eca)
    
}

DERIVATIVE states {
    rates()
    e_q' = rate_e_q 
    f_q' = rate_f_q 
    
}

PROCEDURE rates() {
    LOCAL caConc
    
    caConc = cai
    
    e_timeCourse_t = e_timeCourse_tau ? evaluable
    e_steadyState_x = e_steadyState_rate  / (1 + exp(0 - (v -  e_steadyState_midpoint )/ e_steadyState_scale )) ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=e type=gateHHtauInf), from q10Settings; null
    ? Path not present in component, using factor: 1
    
    e_rateScale = 1 
    
    e_fcond = e_q ^ e_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=e type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)
    e_inf = e_steadyState_x ? path based
    
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=e type=gateHHtauInf), from timeCourse; Component(id=null type=fixedTimeCourse)
    e_tauUnscaled = e_timeCourse_t ? path based
    
    e_tau = e_tauUnscaled  /  e_rateScale ? evaluable
    f_timeCourse_t = f_timeCourse_tau ? evaluable
    f_steadyState_x = f_steadyState_rate  / (1 + exp(0 - (v -  f_steadyState_midpoint )/ f_steadyState_scale )) ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=f type=gateHHtauInf), from q10Settings; null
    ? Path not present in component, using factor: 1
    
    f_rateScale = 1 
    
    f_fcond = f_q ^ f_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=f type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)
    f_inf = f_steadyState_x ? path based
    
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=f type=gateHHtauInf), from timeCourse; Component(id=null type=fixedTimeCourse)
    f_tauUnscaled = f_timeCourse_t ? path based
    
    f_tau = f_tauUnscaled  /  f_rateScale ? evaluable
    h_rateScale = 1 ? evaluable
    h_inf = 1 / (1 + (exp( ( h_ca_half  - caConc) /  h_k ))) ? evaluable
    h_tau = 0 *  h_SEC ? evaluable
    h_q = h_inf ? evaluable
    h_fcond = 1 +(( h_q -1) *  h_alpha ) ? evaluable
    
     
    rate_e_q = ( e_inf  -  e_q ) /  e_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    rate_f_q = ( f_inf  -  f_q ) /  f_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    
}

