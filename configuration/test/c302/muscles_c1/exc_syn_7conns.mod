TITLE Mod file for component: Component(id=exc_syn_7conns type=gradedSynapse)

COMMENT

    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.4.4
         org.neuroml.model   v1.4.4
         jLEMS               v0.9.8.4

ENDCOMMENT

NEURON {
    POINT_PROCESS exc_syn_7conns
    RANGE conductance                       : parameter
    RANGE delta                             : parameter
    RANGE k                                 : parameter
    RANGE Vth                               : parameter
    RANGE erev                              : parameter
    
    RANGE i                                 : exposure
    
    
    NONSPECIFIC_CURRENT i 
    
    RANGE inf                               : exposure
    
    RANGE tau                               : exposure
    POINTER vpeer: derived variable as pointer...
    
    RANGE s_rate                           : conditional derived var
    
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
    
    conductance = 2.3334501E-4 (uS)
    delta = 5 (mV)
    k = 0.025 (kHz)
    Vth = 0 (mV)
    erev = 0 (mV)
}

ASSIGNED {
    ? Standard Assigned variables with baseSynapse
    v (mV)
    celsius (degC)
    temperature (K)
    
    vpeer (mV)                             : derived variable
    
    inf                                    : derived variable
    
    tau (ms)                               : derived variable
    
    i : no units???                        : derived variable
    
    s_rate (kHz)                           : conditional derived var...
    rate_s (/ms)
    
}

STATE {
    s  
    
}

INITIAL {
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    if ((1-  inf  ) <= 1e-4) {
        s = inf ? standard OnCondition
    }
    
    
}

DERIVATIVE states {
    rates()
    s' = rate_s 
    
}

PROCEDURE rates() {
    
    ? DerivedVariable is based on path: peer/v, on: Component(id=exc_syn_7conns type=gradedSynapse), from peer; null
    ? Derived variable: vpeer; its value will be set by a pointer...
    
    inf = 1/(1 + exp((  Vth   - vpeer)/  delta  )) ? evaluable
    tau = (1-  inf  )/  k ? evaluable
    i = -1 * conductance  *  s  * (  erev  -v) ? evaluable
    if ((1-  inf  ) > 1e-4)  { 
        s_rate = (  inf   -   s  )/  tau ? evaluable cdv
    } else  { 
        s_rate = 0 ? evaluable cdv
    }
    
    rate_s = s_rate ? Note units of all quantities used here need to be consistent!
    
     
    
}

