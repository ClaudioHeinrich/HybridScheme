function x = TriInt0(j1,a) 

%Returns the integral of \|x\|^{2a} over the triangle {(x,y) : 0 < x < 0.5, 0 < y < x}
%See Appendix B for details

x=sqrt(2)*j1^(2*a+2)*hypergeom([0.5,1.5+a],1.5,0.5)/(4*(a+1));

end




