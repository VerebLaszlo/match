#include "confuse-parser.h"

int parse_confuse(char const *filename,Parameters *out) {
    int i;
    cfg_t *cfg;
    cfg_t *tmp;

    cfg_opt_t params[] = {
        CFG_FLOAT("m1", 1.0, CFGF_NONE),
        CFG_FLOAT("m2", 1.0, CFGF_NONE),
        CFG_FLOAT("spin1x", 0.1, CFGF_NONE),
        CFG_FLOAT("spin1y", 0.1, CFGF_NONE),
        CFG_FLOAT("spin1z", 0.9, CFGF_NONE),
        CFG_FLOAT("spin2x", 0.11, CFGF_NONE),
        CFG_FLOAT("spin2y", 0.15, CFGF_NONE),
        CFG_FLOAT("spin2z", 0.9, CFGF_NONE),
        CFG_FLOAT("qm1",1, CFGF_NONE),
        CFG_FLOAT("qm2",1, CFGF_NONE),
        CFG_FLOAT("incl",1.4, CFGF_NONE),
        // mi a megfelelője?
        CFG_FLOAT("fI", 50, CFGF_NONE),
        CFG_FLOAT("fL", 2000, CFGF_NONE),
        CFG_FLOAT("fF", 20000, CFGF_NONE),
        CFG_FLOAT("d",1,CFGF_NONE),
        CFG_FLOAT("dt",6E-5,CFGF_NONE),
        CFG_STR("pnorder","twoPointFivePN",CFGF_NONE),
        CFG_STR("spin","QM",CFGF_NONE),
        CFG_END()
    };

    cfg_opt_t root [] = {
        CFG_SEC("params",params,CFGF_MULTI),
        CFG_END()
    };

    cfg = cfg_init(root,CFGF_NONE);
    switch (cfg_parse(cfg,filename)) {
        case CFG_FILE_ERROR:
            printf("Parsing error: %s\n",filename);
            return 1;
            break;
        case CFG_PARSE_ERROR:
            cfg_error(cfg,"Parsing configuration file failed.\n");
            return 1;
            break;
        case CFG_SUCCESS:
            break;
    }

    out->count = cfg_size(cfg,"params");
    out->injParams = calloc(out->count,sizeof(SimInspiralTable));
    out->ppnParams = calloc(out->count,sizeof(PPNParamStruc));


    for (i=0; i < out->count; i++) {
        tmp = cfg_getnsec(cfg,"params",i);
        out->injParams[i].mass1 = cfg_getfloat(tmp,"m1");
        out->injParams[i].mass2 = cfg_getfloat(tmp,"m2");
        out->injParams[i].spin1x = cfg_getfloat(tmp,"spin1x");
        out->injParams[i].spin1y = cfg_getfloat(tmp,"spin1y");
        out->injParams[i].spin1z = cfg_getfloat(tmp,"spin1z");
        out->injParams[i].spin2x = cfg_getfloat(tmp,"spin2x");
        out->injParams[i].spin2y = cfg_getfloat(tmp,"spin2y");
        out->injParams[i].spin2z = cfg_getfloat(tmp,"spin2z");
        out->injParams[i].qmParameter1 = cfg_getfloat(tmp,"qm1");
        out->injParams[i].qmParameter2 = cfg_getfloat(tmp,"qm2");
        out->injParams[i].inclination = cfg_getfloat(tmp,"incl");
        out->injParams[i].f_lower = cfg_getfloat(tmp,"fL");
        out->injParams[i].f_final = cfg_getfloat(tmp,"fF");
        out->injParams[i].distance = cfg_getfloat(tmp,"d");
        out->ppnParams[i].deltaT = cfg_getfloat(tmp,"dt");
    }

    return 0;
}
