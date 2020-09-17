% only yeast
sc = readCbModel(['iMM904.mat']);
sc=setMedia(sc,'C:\Users\user\Desktop\CDM35','CDM35');
sc=changeRxnBounds(sc,'EX_o2_e',0,'l');
% glucose 
glucose_position=getPosOfElementsInArray({'EX_glc__D_e'},sc.rxns);
explicitrxn=getRxn_cobraFormat(sc,glucose_position);
explicitrxn
sc.lb(glucose_position)
% oxygen
oxygen_position=getPosOfElementsInArray({'EX_o2_e'},sc.rxns);
explicitrxn=getRxn_cobraFormat(sc,oxygen_position);
explicitrxn
sc.lb(oxygen_position)
FBAsolution_sc_0 = optimizeCbModel(sc,'max'); % 0

sc.mets = regexprep(sc.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    sc.(fieldToBeCellStr{j})(cellfun(@isempty, sc.(fieldToBeCellStr{j}))) = {''};
end

nameTagsModel = {'SC'};
sc_model = createMultipleSpeciesModel({sc}, nameTagsModel);
sc_model.csense = char('E' * ones(1,numel(sc_model.mets)))';

oxygen_position=getPosOfElementsInArray({'EX_o2[u]'},sc_model.rxns);
sc_model.lb(oxygen_position)=0;

[sc_model.infoCom, sc_model.indCom] = getMultiSpeciesModelId(sc_model, nameTagsModel);
disp(sc_model.infoCom);

% in replacement of YDCom.ub(YDCom.indCom.EXsp(:)) = 1000;
% sc_model.ub(sc_model.indCom.EXsp(find(sc_model.indCom.EXsp(:,2)),2)) = 1000;%
sc_model.ub(sc_model.indCom.EXsp(find(sc_model.indCom.EXsp(:,1)),1)) = 1000;%

sc_model.lb(sc_model.indCom.EXcom)
sc_model.ub(sc_model.indCom.EXcom)

rxnBiomass = strcat(nameTagsModel, sc.rxns(find(sc.c==1)));  % biomass reaction names
rxnBiomassId = findRxnIDs(sc_model, rxnBiomass);  % ids
sc_model.infoCom.spBm = rxnBiomass;  % .spBm for organism biomass reactions
sc_model.indCom.spBm = rxnBiomassId;

printUptakeBoundCom(sc_model, 1);


sc_model.c(rxnBiomassId) =1; 
fbaCom = optimizeCbModel(sc_model);
fbaCom.x(rxnBiomassId) 


options = struct();
options.GRguess = 13;  % initial guess for max. growth rate
options.GRtol = 1e-6;  % tolerance for final growth rate
options.algorithm = 1;  % use the default algorithm (simple guessing for bounds, followed by matlab fzero)



%to make it work for st, I changed line 56 in getCobraSolverParams: 
%Before, it was valDef.feasTol = 1e-9;. Now it is valDef.feasTol = 1e-6;
%basically, this parameter set an upper bound for the largest bound
%violation in the LP problem. 1e-6 is still good. 
% http://www-eio.upc.es/lceio/manuals/cplex-11/html/usrcplex/solveLP20.html

[sol, result] = SteadyCom(sc_model, options);

if ~strcmp(result.stat, 'infeasible')
    for jSp = 1:1
        fprintf('X_%s:  %.6f\n', sc_model.infoCom.spAbbr{jSp}, result.BM(jSp));
    end
    disp(result);
end


% lactococcus only

ll = readCbModel(['iNF517.mat']);

ll=changeRxnBounds(ll,'EX_o2_e',0,'l');
% allow supplement of the media to be adenine,inosine,uracil
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

% if each of these are supplied what is the growth rate?

for i=1:length(missing)
    ll=changeRxnBounds(ll,missing{i},-10,'l');
end
FBAsolution_ll_2 = optimizeCbModel(ll,'max'); % 0.9


ll.mets = regexprep(ll.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ll.(fieldToBeCellStr{j})(cellfun(@isempty, ll.(fieldToBeCellStr{j}))) = {''};
end

nameTagsModel = {'LL'};
ll_model = createMultipleSpeciesModel({ll}, nameTagsModel);
ll_model.csense = char('E' * ones(1,numel(ll_model.mets)))';

[ll_model.infoCom, ll_model.indCom] = getMultiSpeciesModelId(ll_model, nameTagsModel);
disp(ll_model.infoCom);

% in replacement of YDCom.ub(YDCom.indCom.EXsp(:)) = 1000;
% sc_model.ub(sc_model.indCom.EXsp(find(sc_model.indCom.EXsp(:,2)),2)) = 1000;%
ll_model.ub(ll_model.indCom.EXsp(find(ll_model.indCom.EXsp(:,1)),1)) = 1000;%

ll_model.lb(ll_model.indCom.EXcom)
ll_model.ub(ll_model.indCom.EXcom)

rxnBiomass = strcat(nameTagsModel, ll.rxns(find(ll.c==1)));  % biomass reaction names
rxnBiomassId = findRxnIDs(ll_model, rxnBiomass);  % ids
ll_model.infoCom.spBm = rxnBiomass;  % .spBm for organism biomass reactions
ll_model.indCom.spBm = rxnBiomassId;

printUptakeBoundCom(ll_model, 1);


ll_model.c(rxnBiomassId) =1; 
fbaCom = optimizeCbModel(ll_model);
fbaCom.x(rxnBiomassId) 


options = struct();
options.GRguess = 13;  % initial guess for max. growth rate
options.GRtol = 1e-6;  % tolerance for final growth rate
options.algorithm = 1;  % use the default algorithm (simple guessing for bounds, followed by matlab fzero)



%to make it work for st, I changed line 56 in getCobraSolverParams: 
%Before, it was valDef.feasTol = 1e-9;. Now it is valDef.feasTol = 1e-6;
%basically, this parameter set an upper bound for the largest bound
%violation in the LP problem. 1e-6 is still good. 
% http://www-eio.upc.es/lceio/manuals/cplex-11/html/usrcplex/solveLP20.html

[sol, result] = SteadyCom(ll_model, options);

if ~strcmp(result.stat, 'infeasible')
    for jSp = 1:1
        fprintf('X_%s:  %.6f\n', ll_model.infoCom.spAbbr{jSp}, result.BM(jSp));
    end
    disp(result);
end

% Joint Model

