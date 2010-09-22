#include "confuse-parser.h"

void FillFromConfuse(cfg_t *tmp, SimInspiralTable *injParams, PPNParamStruc *ppnParams) {
        injParams->mass1 = cfg_getfloat(tmp,"m1");
        injParams->mass2 = cfg_getfloat(tmp,"m2");
        injParams->spin1x = cfg_getfloat(tmp,"spin1x");
        injParams->spin1y = cfg_getfloat(tmp,"spin1y");
        injParams->spin1z = cfg_getfloat(tmp,"spin1z");
        injParams->spin2x = cfg_getfloat(tmp,"spin2x");
        injParams->spin2y = cfg_getfloat(tmp,"spin2y");
        injParams->spin2z = cfg_getfloat(tmp,"spin2z");
        injParams->qmParameter1 = cfg_getfloat(tmp,"qm1");
        injParams->qmParameter2 = cfg_getfloat(tmp,"qm2");
        injParams->inclination = cfg_getfloat(tmp,"incl");
        injParams->f_lower = cfg_getfloat(tmp,"fL");
        injParams->f_final = cfg_getfloat(tmp,"fF");
        injParams->distance = cfg_getfloat(tmp,"d");
        ppnParams->deltaT = cfg_getfloat(tmp,"dt");
}

int parse_confuse(char const *filename,Parameters *out) {
    int i;
    cfg_t *cfg, *tmp;
    SimInspiralTable default_siminsp;
    PPNParamStruc default_ppn;

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
        // mi a megfelelője?
        CFG_FLOAT("d",1,CFGF_NONE),
        CFG_FLOAT("dt",6E-5,CFGF_NONE),
        CFG_STR("pnorder","twoPointFivePN",CFGF_NONE),
        CFG_STR("spin","QM",CFGF_NONE),
        CFG_END()
    };

    cfg_opt_t root [] = {
        CFG_SEC("waveDefault",params,CFGF_NONE),
        CFG_SEC("waveParams",params,CFGF_MULTI | CFGF_TITLE),
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

    tmp = cfg_getsec(cfg,"waveDefault");
    FillFromConfuse(tmp,&default_siminsp, &default_ppn);
    cfg_free(cfg);

    cfg_opt_t params_def[] = {
        CFG_FLOAT("m1", default_siminsp.mass1, CFGF_NONE),
        CFG_FLOAT("m2", default_siminsp.mass2, CFGF_NONE),
        CFG_FLOAT("spin1x", default_siminsp.spin1x, CFGF_NONE),
        CFG_FLOAT("spin1y", default_siminsp.spin1y, CFGF_NONE),
        CFG_FLOAT("spin1z", default_siminsp.spin1z, CFGF_NONE),
        CFG_FLOAT("spin2x", default_siminsp.spin2x, CFGF_NONE),
        CFG_FLOAT("spin2y", default_siminsp.spin2y, CFGF_NONE),
        CFG_FLOAT("spin2z", default_siminsp.spin2z, CFGF_NONE),
        CFG_FLOAT("qm1",default_siminsp.qmParameter1, CFGF_NONE),
        CFG_FLOAT("qm2",default_siminsp.qmParameter2, CFGF_NONE),
        CFG_FLOAT("incl",default_siminsp.inclination, CFGF_NONE),
        CFG_FLOAT("fI", 50, CFGF_NONE),
        CFG_FLOAT("fL", default_siminsp.f_lower, CFGF_NONE),
        CFG_FLOAT("fF", default_siminsp.f_final, CFGF_NONE),
        CFG_FLOAT("d",1,CFGF_NONE),
        CFG_FLOAT("dt",default_ppn.deltaT,CFGF_NONE),
        CFG_STR("pnorder","twoPointFivePN",CFGF_NONE),
        CFG_STR("spin","QM",CFGF_NONE),
        CFG_END()
    };

    cfg_opt_t root_def [] = {
        CFG_SEC("waveDefault",params,CFGF_NONE),
        CFG_SEC("waveParams",params_def,CFGF_MULTI | CFGF_TITLE),
        CFG_END()
    };

    cfg = cfg_init(root_def,CFGF_NONE);
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

    out->count = cfg_size(cfg,"waveParams");
    out->title = calloc(out->count,sizeof(*out->title));
    out->injParams = calloc(out->count,sizeof(SimInspiralTable));
    out->ppnParams = calloc(out->count,sizeof(PPNParamStruc));


    for (i=0; i < out->count; i++) {
        tmp = cfg_getnsec(cfg,"waveParams",i);
        FillFromConfuse(tmp,&(out->injParams[i]),&(out->ppnParams[i]));
        out->title[i] = strdup(cfg_title(tmp));
    }

    cfg_free(cfg);

    return 0;
}
