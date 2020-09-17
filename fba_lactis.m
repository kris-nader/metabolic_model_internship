ll = readCbModel(['iNF517.mat']);

%% restrict o2
ll=changeRxnBounds(ll,'EX_o2_e',0,'l');

%% find missing nutrients
FBAsolution_ll_0 = optimizeCbModel(ll,'max'); %0.0426
temp=findExcRxns(ll);
excRxns=ll.rxns(temp);
negative_flux=ll.rxns(find(FBAsolution_ll_0.x<0));
nutrients=intersect(excRxns,negative_flux);
positions=getPosOfElementsInArray(nutrients,ll.rxns);
lb_nutrients=ll.lb(positions);

ll=setMedia(ll,'C:\Users\user\Desktop\CDM35','CDM35');
FBAsolution_ll_1 = optimizeCbModel(ll,'max'); %0
allowedbackwards=ll.rxns(find(ll.lb<0));
temp=findExcRxns(ll);
excRxns=ll.rxns(temp);
allowedNutrients=intersect(allowedbackwards,excRxns);
missing=setdiff(nutrients,allowedNutrients);

