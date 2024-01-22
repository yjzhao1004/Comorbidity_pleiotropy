load("AHBA_Mean_reanote.mat")
genes=probeInformation.EntrezID; 
geneSymbol = probeInformation.GeneSymbol;

unique_gene = readtable("unique_pleiogeneid.txt");
unique_gene = table2array(unique_gene);

[tf, index] = ismember(unique_gene,genes);
index(find(index==0))=[];

cort_expMeanRaw = cell(6,1);
for subject=1:6
     if subject ==1
         exp_raw = expressionMean{subject}(:,index); 
     elseif subject ==2
         exp_raw = expressionMean{subject}(:,index);
     elseif subject ==3 
         exp_raw = expressionMean{subject}(:,index);
     elseif subject ==4
         exp_raw = expressionMean{subject}(:,index);
     elseif subject ==5
         exp_raw = expressionMean{subject}(:,index);
     elseif subject ==6 
         exp_raw = expressionMean{subject}(:,index);
     end
     cort_expMeanRaw{subject} = exp_raw;
end

%scaledRobustSigmoid
cort_expMeanScaled = cell(6,1);
cort_badNorm = cell(6,1);%badly-normalised genes

%first normalized with sample across gene, then normalized across sample
for subject = 1:6
    exp_raw = cort_expMeanRaw{subject};
    exp_norm1 = BF_NormalizeMatrix(exp_raw','scaledRobustSigmoid')';
    exp_norm2 = BF_NormalizeMatrix(exp_norm1,'scaledRobustSigmoid');
    temp  = find(any(isnan(exp_norm2),1));
    cort_expMeanScaled{subject} = exp_norm2;
    cort_badNorm{subject} = temp;
end


badNorm = horzcat(cort_badNorm{1},cort_badNorm{2},cort_badNorm{3},...
    cort_badNorm{4},cort_badNorm{5},cort_badNorm{6});
ubadNorm = unique(badNorm);

cort_l_expMS = vertcat(cort_expMeanScaled{1},cort_expMeanScaled{2},cort_expMeanScaled{3},...
    cort_expMeanScaled{4},cort_expMeanScaled{5},cort_expMeanScaled{6});


%remove the bad normalized genes
cort_l_expMS(:,ubadNorm) = [];
index(ubadNorm) = [];



gene_symbol_selected = geneSymbol(index);
gene_id_selected = genes(index);



save AHBA_Mean_scaled_reanote_selectedpleiogenes cort_expMeanScaled cort_l_expMS...
    gene_id_selected gene_symbol_selected Coordinatesall;


load("AHBA_Mean_scaled_reanote_selectedpleiogenes.mat")
cort_l_expMS_meangenes = mean(cort_l_expMS,2);


parcelExpressionselectedgenes = parcelExpression(:,index);
parcelExpressionselectedgenes_meangenes = mean(parcelExpressionselectedgenes,2);

parcelExpressionselectedgenes_meangenes = parcelExpression(:,find(gene_symbol_selected=="PITPNM2"));

parcelExpressionselectedgenes_meangenes_plot = [parcelExpressionselectedgenes_meangenes',zeros(1,34)];

corvol_life{:,1}=parcelExpressionselectedgenes_meangenes_plot;

volname = CBIG_text2cell('vol_name.txt');
corvol_life{:,2} = volname(1:68);


covv_life_d=corvol_life{:,1};
covv_life_d=array2table(covv_life_d);

covv_life_d.Properties.VariableNames=corvol_life{:,2}';

min(parcelExpressionselectedgenes_meangenes_plot)
max(parcelExpressionselectedgenes_meangenes_plot)

CV1_d_fsa5 = parcel_to_surface(parcelExpressionselectedgenes_meangenes_plot, 'aparc_fsa5');
f = figure,
plot_cortical(CV1_d_fsa5, 'surface_name', 'fsa5', 'color_range', ...
[0.25 0.65], 'cmap', 'Reds')





