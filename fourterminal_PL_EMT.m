%N*M resistors,each with 8 variables:i1,i2,i3,i4,j1,j2,j3,j4
%N*M*8 variables
%N*M*8 equations

function[r_eq]=fourterminal_PL_EMT(N,M,nsd,nmean,H)

%constants

    Vx=1;

    rng(0,'twister')
    musd = 0;
    mumean = 1;
    
    muValues=musd.*randn(1,N*M) + mumean;
        
    nValues = nsd.*randn(1,N*M) + nmean;
    
    rhoValues = abs(1./(nValues.*muValues));
    
    t=1;
    rValues=rhoValues./(pi*t);
    BValues = muValues.*H;
    
    phi=0.14;
    g=1/phi;
    
    %a=-g+(pi/4)*B;
        aValues= -g+(pi/4).*BValues;
    %b=g+(pi/4)*B;
        bValues= g+(pi/4).*BValues;
    %c=0.35-(pi/4)*B;
        cValues=0.35-(pi/4).*BValues;
    %d=-0.35-(pi/4)*B;
        dValues=-0.35-(pi/4).*BValues;

%coefficients matrix
    CoeffsMatrix=zeros(8*N*M,8*N*M);

%constants matrix
    ConstantsMatrix=zeros(8*N*M,1);
    ConstantsMatrix(1+2*N*(M-1)+2*M*(N-1)+2*M+N:2*N*(M-1)+2*M*(N-1)+2*M+2*N)=Vx;

%current matrix
    %i1Coeffs=CoeffsMatrix(:,1:N*M);
    %i2Coeffs=CoeffsMatrix(:,N*M+1:2*N*M);
    %i3Coeffs=CoeffsMatrix(:,2*N*M+1:3*N*M);
    %i4Coeffs=CoeffsMatrix(:,3*N*M+1:4*N*M);

%voltage matrix
    %v1Coeffs=CoeffsMatrix(:,4*N*M+1:5*N*M)
    %v2Coeffs=CoeffsMatrix(:,5*N*M+1:6*N*M)
    %v3Coeffs=CoeffsMatrix(:,6*N*M+1:7*N*M)
    %v1Coeffs=CoeffsMatrix(:,7*N*M+1:8*N*M)
    
%equations

%x direction equations
    %ix equations: i3 of resistor = -i1 of next resistor in each row
    %There are N*(M-1) ix equations

        m=2*N*M+1;
        for j=1:N
            for n=(j-1)*(M-1)+1:j*(M-1)
                    CoeffsMatrix(n,m)=1;
                    CoeffsMatrix(n,m-2*N*M+1)=1;
              m=m+1;
            end
            m=m+1;
        end
            
    %vx equations: v3 of resistor = v1 of next resistor in each row
    %There are N*(M-1) vx equations

        m=6*N*M+1;
        for j=1:N
            for n=1+(j-1)*(M-1)+N*(M-1):j*(M-1)+N*(M-1)
                    CoeffsMatrix(n,m)=1;
                    CoeffsMatrix(n,m-2*N*M+1)=-1;
              m=m+1;
            end
            m=m+1;
        end
            
%y direction equations
    %iy equations: i4 of resistor=-i2 of next resistor in each column
    %There are M*(N-1) iy equations
    
         m=3*N*M+1;
            for n=2*N*(M-1)+1:2*N*(M-1)+M*(N-1)
                CoeffsMatrix(n,m)=1;
                CoeffsMatrix(n,m-2*N*M+M)=1;
                m=m+1;
            end
            
    %vy equations: v4 of resistor=v2 of next resistor in each column
    %There are M*(N-1) iy equations
    
         m=7*N*M+1;
            for n=2*N*(M-1)+M*(N-1)+1:2*N*(M-1)+2*M*(N-1)
                CoeffsMatrix(n,m)=1;
                CoeffsMatrix(n,m-2*N*M+M)=-1;
                m=m+1;
            end
            
     %i top equations: i2 of top row = 0
     %there are M itop equations
     
        m=N*M+1;
            for n=(2*N*(M-1)+2*M*(N-1))+1:(2*N*(M-1)+2*M*(N-1)+M)
                CoeffsMatrix(n,m)=1;
                m=m+1;
            end

     %i bottom equations: i4 of bottom row = 0
     %there are M ibottom equations
        
        m=4*N*M-M+1;
            for n=2*N*(M-1)+2*M*(N-1)+M+1:2*N*(M-1)+2*M*(N-1)+2*M
                CoeffsMatrix(n,m)=1;
                m=m+1;
            end
            
%V equations

    %V left equation: v1 on left edge = 0
    %There are N VLeftEquations

        m=4*N*M+1;
        for n=1+4*N*M-2*N:4*N*M-N
            CoeffsMatrix(n,m)=1;
            m=m+M;
        end
        
    %V right equation: v3 on right edge = Vx
    %There are N VRightEquations
    
        m=6*N*M+M;
        for n=1+2*N*(M-1)+2*M*(N-1)+2*M+N:2*N*(M-1)+2*M*(N-1)+2*M+2*N
            CoeffsMatrix(n,m)=1;
            m=m+M;
        end
        
%KCL equations: i1+i2+i3+i4 at each resistor=0
    %there are NM KCL equations
    
    m=1;
    for n=4*M*N+1:5*M*N
        CoeffsMatrix(n,m)=1; %i1
        CoeffsMatrix(n,m+N*M)=1; %i2
        CoeffsMatrix(n,m+2*N*M)=1; %i3
        CoeffsMatrix(n,m+3*N*M)=1; %i4
        m=m+1;
    end
        
        
%Resistance Matrix equations
%     %There are 4*N*M RMatrix equations 

            m=1;
            for n=5*M*N+1:3:8*M*N
                CoeffsMatrix(n,m)=rValues(m)*aValues(m); %first i1
                CoeffsMatrix(n,m+N*M)=rValues(m)*bValues(m); %first i2
                CoeffsMatrix(n,m+2*N*M)=rValues(m)*cValues(m); %first i3
                CoeffsMatrix(n,m+3*N*M)=rValues(m)*dValues(m); %first i4
                CoeffsMatrix(n,m+4*N*M)=1; %v1
                CoeffsMatrix(n,m+5*N*M)=-1; %v2
                
                CoeffsMatrix(n+1,m)=rValues(m)*dValues(m); %second i1
                CoeffsMatrix(n+1,m+N*M)=rValues(m)*aValues(m); %second i2
                CoeffsMatrix(n+1,m+2*N*M)=rValues(m)*bValues(m); %second i3
                CoeffsMatrix(n+1,m+3*N*M)=rValues(m)*cValues(m); %second i4
                CoeffsMatrix(n+1,m+5*N*M)=1; %v2
                CoeffsMatrix(n+1,m+6*N*M)=-1; %v3
                
                CoeffsMatrix(n+2,m)=rValues(m)*cValues(m); %third i1
                CoeffsMatrix(n+2,m+N*M)=rValues(m)*dValues(m); %third i2
                CoeffsMatrix(n+2,m+2*N*M)=rValues(m)*aValues(m); %third i3
                CoeffsMatrix(n+2,m+3*N*M)=rValues(m)*bValues(m); %third i4
                CoeffsMatrix(n+2,m+6*N*M)=1; %v3
                CoeffsMatrix(n+2,m+7*N*M)=-1; %v4
                
%                 CoeffsMatrix(n+3,m)=rValues(m)*bValues(m); %fourth i1
%                 CoeffsMatrix(n+3,m+N*M)=rValues(m)*cValues(m); %fourth i2
%                 CoeffsMatrix(n+3,m+2*N*M)=rValues(m)*dValues(m); %fourth i3
%                 CoeffsMatrix(n+3,m+3*N*M)=rValues(m)*aValues(m); %fourth i4
%                 CoeffsMatrix(n+3,m+7*N*M)=-1; %v4
%                 CoeffsMatrix(n+3,m+4*N*M)=1; %v1
                m=m+1;
            end
                
%FINAL CALCULATIONS

S=sparse(CoeffsMatrix);
Values=S\ConstantsMatrix;
% Values=CoeffsMatrix\ConstantsMatrix;
SumInputI=sum(Values(2*N*M+M:M:3*N*M));
r_eq=Vx/SumInputI;

% Values

end     
