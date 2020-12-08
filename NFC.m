function [NFC_values] = NFC(ORF , RC_profiles)
% calc the NFC value for each gene in ORF relative to RC profiles
% NFC is calculated by read count (in RC profiles) divided by the average
% read count for the gene
NFC_values = {};
geneNFC = [];
for i = [1 : length(ORF)]
    % for every codon in every gene
    geneNTseq = char(ORF(i));
    geneRCprofile = cell2mat(RC_profiles (i));
    NFC_AVG = mean (geneRCprofile);
    for j = [1:3:length(geneNTseq)]
%         codon = geneNTseq(j:j+2);
        codonRC = geneRCprofile ((j+2)/3);
        codonNFC = codonRC/NFC_AVG;
        geneNFC = horzcat(geneNFC,codonNFC);
    end 
    NFC_values = horzcat (NFC_values,geneNFC);
end
    
end

