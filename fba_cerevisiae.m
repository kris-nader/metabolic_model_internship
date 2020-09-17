
sc = readCbModel(['iMM904.mat']);
% no cdm35 +o2+ glucose
sc_model1=sc;
sc_model1 = changeRxnBounds(sc_model1,'EX_glc__D_e',-10,'l');
sc_model1  =changeRxnBounds(sc_model1,'EX_o2_e',-10,'l');
FBAsolution_sc_1 = optimizeCbModel(sc_model1,'max'); %0.5365

% cdm35 + glucose +o2
sc_model2=sc;
sc_model2=setMedia(sc_model2,'C:\Users\user\Desktop\CDM35','CDM35');
sc_model2 = changeRxnBounds(sc_model2,'EX_glc__D_e',-10,'l');
sc_model2  =changeRxnBounds(sc_model2,'EX_o2_e',-10,'l');
FBAsolution_sc_2 = optimizeCbModel(sc_model2,'max'); %0.7074

% cdm35 + glucose -o2
sc_model3=sc;
sc_model3=setMedia(sc_model3,'C:\Users\user\Desktop\CDM35','CDM35');
sc_model3 = changeRxnBounds(sc_model3,'EX_glc__D_e',-10,'l');
sc_model3  =changeRxnBounds(sc_model3,'EX_o2_e',0,'l');
FBAsolution_sc_3 = optimizeCbModel(sc_model3,'max'); %0

% cdm35- glucose +o2
sc_model4=sc;
sc_model4=setMedia(sc_model4,'C:\Users\user\Desktop\CDM35','CDM35');
sc_model4 = changeRxnBounds(sc_model4,'EX_glc__D_e',0,'l');
sc_model4  =changeRxnBounds(sc_model4,'EX_o2_e',-10,'l');
FBAsolution_sc_4 = optimizeCbModel(sc_model4,'max'); %0.4017

% cdm35 -glucose -o2
sc_model5=sc;
sc_model5=setMedia(sc_model5,'C:\Users\user\Desktop\CDM35','CDM35');
sc_model5 = changeRxnBounds(sc_model5,'EX_glc__D_e',0,'l');
sc_model5  =changeRxnBounds(sc_model5,'EX_o2_e',0,'l');
FBAsolution_sc_5 = optimizeCbModel(sc_model5,'max'); %0

% no cdm35 +glucose -o2
sc_model7=sc;

sc_model7  =changeRxnBounds(sc_model7,'EX_o2_e',0,'l');
FBAsolution_sc_8 = optimizeCbModel(sc_model7,'max'); %0



%% which nutrients can be provided by yeast
sc_model6=sc;

fba_temp=[];
for i=1:length(missing)
    sc_model6=changeRxnBounds(sc,missing{i},0.01,'l');
    fba_temp(i)=optimizeCbModel(sc_model6,'max').f;
    sc_model6=sc;
end

%%
% Uptake reaction
model=changeRxnBounds(model,'EX_glc__D_e',-10,'l');

% Force the model to produce and export
model=changeRxnBounds(model,'EX_glc__D_e',+10,'l');

% Prevent uptake
model=changeRxnBounds(model,'EX_glc__D_e',0,'l');


sc = readCbModel(['iMM904.mat']);
% no cdm35 +o2+ glucose
sc_model1=sc;
sc_model1 = changeRxnBounds(sc_model1,'EX_glc__D_e',-10,'l');
sc_model1  =changeRxnBounds(sc_model1,'EX_o2_e',0,'l');
FBAsolution_sc_1 = optimizeCbModel(sc_model1,'max'); %0.0



