function [NFCmap] = codon_NFC_map(ORF, NFCval)
% mapping NFC for each codon
% initializing struct
codon_names = fieldnames(codoncount(char(ORF(1))));
NFCmap = ([]);

for i = [1:length(codon_names)]
    NFCmap.(char(codon_names(i))) = [];
end

for i = 1 : length(ORF)
    % for every codon in every gene
    geneNTseq = char(ORF(i));
    gene_NFC = NFCval{1,i};
    for j = 1:3:length(geneNTseq)
        codonNT = geneNTseq(j:j+2);
        codonNFC = gene_NFC((j+2)/3);
        NFCmap.(codonNT) = horzcat(NFCmap.(codonNT),codonNFC);
    end 

end

for i = 1:length(codon_names)
    codon = char(codon_names(i));
    NFC_values_vect = NFCmap.(codon);
    max_NFC = max(NFC_values_vect);
    norm_NFC_val = NFC_values_vect./max_NFC;
    NFCmap.(codon) = norm_NFC_val;
end
    

end

