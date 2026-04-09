function [Th_out, Tc_out, epsilon] = CounterFlowFunction(mhdot, mcdot, ch, cc, Thi, Tci, U, L, D)
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


N=10000;
dx=L/N;

Tc=zeros(1,N+1);
Th=zeros(1,N+1);
Tc(N+1)=Tci;
Th(1)=Thi;

if (mhdot*ch==mcdot*cc) %%in this case Tlm=DT1=DT2
        Tco=(Tci+U*pi*D*L*Thi/(mcdot*cc))/(1+U*pi*D*L/(mcdot*cc));
        Tho=Thi+mcdot*cc*(Tci-Tco)/(mhdot*ch);

        Tc(1)=Tco;
        Th(N+1)=Tho;

        for i=2:N
            Tci=Tco+U*pi*D*(i*dx)*(Tco-Thi)/(mcdot*cc);
            Tc(i)=Tci;
            Th(i)=Thi+mcdot*cc*(Tci-Tco)/(mhdot*ch);
        end
        
else
        Tco2=Thi;       %Start by finding Tco and Tho
        Tco1=Tci;       
        F=1;
        Fold=0;
        while abs(F)>0.00000001
            Tco=Tco1+(Tco2-Tco1)/2;
            Tho=Thi+(mcdot*cc)*(Tci-Tco)/(mhdot*ch);
            DT1=Thi-Tco;
            DT2=Tho-Tci;
            Tlm=(DT1-DT2)/log(DT1/DT2);
            F=mcdot*cc*(Tco-Tci)-U*pi*D*L*Tlm;
            if (F == Fold)
                break;
            elseif (F<0)
                Tco1=Tco;    
            else
                Tco2=Tco; 
            end
            
            Fold = F;
        end
    
        Tc(1)=Tco;
        Th(N+1)=Tho;
   
        if (mhdot*ch < mcdot*cc) 
            for i=1:N    %Now march through heat exchanger
                
                Tci2=Tco;   %Each dx of exchanger, we will know the local Thi and Tco
                Tci1=Tc(N+1);   %Iterate to solve for the local Tci, and then compute the local Tho
                
                F=1;
                Fold=0;
                while abs(F)>0.000000001
                    Tci=Tci1+(Tci2-Tci1)/2;
                    Tho=Thi+(mcdot*cc)*(Tci-Tco)/(mhdot*ch);
                    DT1=Thi-Tco;
                    DT2=Tho-Tci;
                    Tlm=(DT1-DT2)/log(DT1/DT2);
                    F=mcdot*cc*(Tco-Tci)-U*pi*D*dx*Tlm;
                    if (F == Fold || isnan(F))
                        break;
                    elseif (F>0)
                        Tci1=Tci;
                    else
                        Tci2=Tci;
                    end
                    Fold=F;
                end

                Th(i+1)=Tho;
                Tc(i+1)=Tci;
                
                Thi=Tho;
                Tco=Tci;
                

            end
            
            
        else  
            for i=N:-1:1    %Now march through heat exchanger
                
                Thi2=Th(1);   %Each dx of exchanger, we will know the local Thi and Tco
                Thi1=Tho;   %Iterate to solve for the local Tci, and then compute the local Tho
                
                F=1;
                Fold=0;
                while abs(F)>0.000000001
                    Thi=Thi1+(Thi2-Thi1)/2;
                    Tco=Tci+(mhdot*ch)*(Thi-Tho)/(mcdot*cc);
                    DT1=Thi-Tco;
                    DT2=Tho-Tci;
                    Tlm=(DT1-DT2)/log(DT1/DT2);
                    F=mhdot*ch*(Thi-Tho)-U*pi*D*dx*Tlm;
                    if (F == Fold || isnan(F))
                        break;
                    elseif (F<0)
                        Thi1=Thi;
                    else
                        Thi2=Thi;
                    end
                    Fold=F;
                end
 
                Th(i)=Thi;
                Tc(i)=Tco;
                
                Tho=Thi;
                Tci=Tco;

            end
            
        end
end

x=0:L/N:L;
plot(x,Th,'r',x,Tc,'b','LineWidth',2);
axis([0 L yl yh]);
set(gca,'FontSize',12);
ylabel(sprintf('Temperature (%cC)',char(176)),'FontSize',14)
xlabel('Length (m)','FontSize',14)

Th_out = Th(end);
Tc_out = Tc(1);

    %Compute effectiveness
    qh=(mhdot*ch)*(Th(1)-Th_out);
    if (mhdot*ch)<(mcdot*cc)
        qmax=(mhdot*ch)*(Th(1)-Tc(N+1));
    else
        qmax=(mcdot*cc)*(Th(1)-Tc(N+1));
    end
    epsilon=qh/qmax;

end

