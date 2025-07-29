function x = inverseFourierTransformVector(X,w,t)

z = exp(1i*w(:));
zz = z.^(-t(:)');
x = zz \ X;


