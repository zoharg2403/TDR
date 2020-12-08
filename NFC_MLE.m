function [MLE,mu,sigma] = NFC_MLE(NFCmap)
codon_names = fieldnames(NFCmap);

MLE = ([]);
mu = ([]);
sigma = ([]);

log_likelihood =0;
for i = 1:length(codon_names)
    codon = char(codon_names (i));
    codonNFC = NFCmap.(codon); 
    codonNFC = codonNFC(codonNFC~=0);
    
    codon_mu = mean(codonNFC);  
    codon_sigma = std(codonNFC);
    codon_x = 0:0.001:3;

    log_dist_func = mle(codonNFC,'pdf',@(codon_x,codon_mu,codon_sigma)log_normal_pdf(codon_x,codon_mu,codon_sigma),'start',[1,1]);
    log_likelihood = log_likelihood + log_dist_func;
    MLE.(codon) = log_likelihood;
    mu.(codon) = log_likelihood(1);
    sigma.(codon) = log_likelihood(2);
    log_likelihood = 0;
end
end

