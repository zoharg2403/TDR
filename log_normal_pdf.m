function  [value] = log_normal_pdf(x, mu, sigma)

value=(1./(x*sigma*sqrt(2*pi))).*exp(-((log(x)-mu).^2)./(2*sigma^2));

value(value <= 0) = eps/2;

end

