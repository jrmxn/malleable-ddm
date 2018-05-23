sr = ddm_def('');
mk = sr.getmodelKey(0);
%%
id_model = debi_model(0,'de','bi');
id_model(mk.s) = 1;
id_model(mk.a) = 1;
id_model(mk.t) = 1;
id_model(mk.v) = 1;
id_model = debi_model(id_model,'bi','de');

id_fit = debi_model(0,'de','bi');
id_fit(mk.s) = 0;
id_fit(mk.a) = 1;
id_fit(mk.t) = 1;
id_fit(mk.v) = 1;
id_fit = debi_model(id_fit,'bi','de');

sr.subject = 'PRPAE';
sr.path_data = 'testing.csv';

sr = sr.ddm_init(id_model,id_fit);
sr = sr.ddm_search;