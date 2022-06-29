function result = output(lambda, I, G1, G2, a,b,c,d,e,g,x,y,z)
    eta = 1.2;
    p=0.05; %Mbar
    result = I*((1-lambda).^b)*((eta-1-a).^x) + G1*((1-lambda).^c)*(lambda.^d)*(p.^y) + G2*((1-lambda).^e)*(lambda.^g)*(p.^z);
end
