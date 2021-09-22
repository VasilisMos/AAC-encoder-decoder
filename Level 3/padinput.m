function xout = padinput(x,PRINTINFO)
%Padds the input at the beginning and at the end with zeros so that total sample are divided perfectly with 2048
Ntot = max(size(x));
temp = ceil(Ntot/2048)*2048;
prepad = 1024;
postpad = temp-Ntot + 1024;

if(PRINTINFO)fprintf("Padding input, Noriginal = %d, Nnew=%d prepad=%d,postpad=%d\n",Ntot,temp,prepad,postpad);end
ch1 = x(:,1);
ch2 = x(:,2);

ch1 = padarray(ch1,prepad,0,"pre");
ch1 = padarray(ch1,postpad,0,"post");

ch2 = padarray(ch2,prepad,0,"pre");
ch2 = padarray(ch2,postpad,0,"post");

xout = [ch1 ch2];

end