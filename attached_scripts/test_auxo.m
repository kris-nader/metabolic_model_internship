
function dispResults = testListOfAuxotrophies
model = readSBML('C:\Users\user\Desktop\saccharomyces_cerevisiae\iMM904.xml',1000);
listOfReactions = model.rxns(find(findExcRxns(model)));
results = zeros(length(listOfReactions),1);
for i = 1:length(listOfReactions)
 result = testAuxotrophy(model, listOfReactions{i});
 results(i) = result;
end
dispResults = {listOfReactions num2cell(results)}
end

function isAbleToGrow = testAuxotrophy(model,rxnID)
%function description
isAbleToGrow = 0;
model = changeRxnBounds(model,rxnID,0,'l');
fba = optimizeCbModel(model);
if fba.f > 0.001
    isAbleToGrow = 1;
end
end


% credit to Jessica Singh
