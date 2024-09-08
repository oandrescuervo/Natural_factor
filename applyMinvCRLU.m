function s=applyMinvCRLU(r,L,U,D)
s1=U'\r;
s2=L'\s1;
s3=D\s2;
s4=L\s3;
s=U\s4;

