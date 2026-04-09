function [Th_out, Tc_out, epsilon] = ParallelFlowFunction(mhdot, mcdot, ch, cc, Thi, Tci, U, L, D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   mhdot - mass flow rate of hot fluid
%   mcdot - mass flow rate of cool fluid
%   ch - specific heat of hot fluid
%   cc - specific heat of cold fluid
%   Thi - temperature at the hot fluid inlet
%   Tci - temperature at the cold fluid inlet
%   U - the overall heat transfer coefficient
%   L - the length of the heat exchanger
%   D - mean diameter of the tube separating the fluids
    
   yh = Thi+(Thi-Tci)/10;
   yl = Tci-(Thi-Tci)/10;

dx=L/10000;

Tc=zeros(1,L/dx+1);
Th=zeros(1,L/dx+1);
Tc(1)=Tci;
Th(1)=Thi;

    for i=1:L/dx
        
        Tco2=Thi;
        Tco1=Tci;
          
       % disp(i)
        
        F=1;
        while abs(F)>0.00001
            Tco=Tco1+(Tco2-Tco1)/2; 
            F=mcdot*cc*(Tco-Tci)-U*pi*D*dx*((Thi-Tci)-(Thi+mcdot*cc*(Tci-Tco)/(mhdot*ch)-Tco))/log((Thi-Tci)/(Thi+mcdot*cc*(Tci-Tco)/(mhdot*ch)-Tco));
            if (F<0)
                Tco1=Tco;    
            else
                Tco2=Tco; 
            end
        end
    
        Thi=Thi+mcdot*cc*(Tci-Tco)/(mhdot*ch);
        Tci=Tco;
        
        Th(i+1)=Thi;
        Tc(i+1)=Tci;
    end
    
    x=0:L/10000:L;
    plot(x,Th,'r',x,Tc,'b','LineWidth',2);
    axis([0 L yl yh]);
    set(gca,'FontSize',12);
    ylabel(sprintf('Temperature (%cC)',char(176)),'FontSize',14)
    xlabel('Length (m)','FontSize',14)
    
    Th_out = Th(end);
    Tc_out = Tc(end);
    
    %Compute effectiveness
    qh=(mhdot*ch)*(Th(1)-Th_out);
    if (mhdot*ch)<(mcdot*cc)
        qmax=(mhdot*ch)*(Th(1)-Tc(1));
    else
        qmax=(mcdot*cc)*(Th(1)-Tc(1));
    end
    epsilon=qh/qmax;
    

end

