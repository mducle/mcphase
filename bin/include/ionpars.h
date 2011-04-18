// declare cfield  module
extern "C" void getpar(char * iontype, int * dj, double * alpha, double * beta, double * gamma, double * lande,
       double * rh2, double * rh4,double * rh6, int *nof_electrons );

extern "C" void cfield_mcphasnew(char * iontype,double *Jxr,double * Jxi,  double * Jyr, double * Jyi, double * Jzr, double * Jzi,
                              double * mo22sr, double * mo22si,
                              double * mo21sr, double * mo21si,
                              double * mo20cr, double * mo20ci,
                              double * mo21cr, double * mo21ci,
                              double * mo22cr, double * mo22ci,

                              double * mo33sr, double * mo33si,
                              double * mo32sr, double * mo32si,
                              double * mo31sr, double * mo31si,
                              double * mo30cr, double * mo30ci,
                              double * mo31cr, double * mo31ci,
                              double * mo32cr, double * mo32ci,
                              double * mo33cr, double * mo33ci,

                              double * mo44sr, double * mo44si,
                              double * mo43sr, double * mo43si,
                              double * mo42sr, double * mo42si,
                              double * mo41sr, double * mo41si,
                              double * mo40cr, double * mo40ci,
                              double * mo41cr, double * mo41ci,
                              double * mo42cr, double * mo42ci,
                              double * mo43cr, double * mo43ci,
                              double * mo44cr, double * mo44ci,

                              double * mo55sr, double * mo55si,
                              double * mo54sr, double * mo54si,
                              double * mo53sr, double * mo53si,
                              double * mo52sr, double * mo52si,
                              double * mo51sr, double * mo51si,
                              double * mo50cr, double * mo50ci,
                              double * mo51cr, double * mo51ci,
                              double * mo52cr, double * mo52ci,
                              double * mo53cr, double * mo53ci,
                              double * mo54cr, double * mo54ci,
                              double * mo55cr, double * mo55ci,

                              double * mo66sr, double * mo66si,
                              double * mo65sr, double * mo65si,
                              double * mo64sr, double * mo64si,
                              double * mo63sr, double * mo63si,
                              double * mo62sr, double * mo62si,
                              double * mo61sr, double * mo61si,
                              double * mo60cr, double * mo60ci,
                              double * mo61cr, double * mo61ci,
                              double * mo62cr, double * mo62ci,
                              double * mo63cr, double * mo63ci,
                              double * mo64cr, double * mo64ci,
                              double * mo65cr, double * mo65ci,
                              double * mo66cr, double * mo66ci,
                              double * modxcr, double * modxci,
                              double * modycr, double * modyci,
                              double * modzcr, double * modzci,int * dimj,
                              double * alpha,double * beta, double * gamma, double * gJ,
                              double * rh2, double * rh4,double * rh6, int * nof_electrons);


