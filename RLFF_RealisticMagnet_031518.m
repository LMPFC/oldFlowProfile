%function RLFF_OmegaOptimize_013118

clear
clc
close all

%% Define constants
pm = 10;                    % Magic magnet number
mu0 = 4*pi*1e-7;            % Permittivity of free space


%% Magnetic field generation
% radialRes = 400;
% lengthRes = 800;
% L = 0.0254/4;
% rhoMaxVal = 0.1;
% zMaxVal = 0.1;
% rhoSpace = linspace(0,rhoMaxVal,radialRes);
% zSpace = linspace(-zMaxVal,zMaxVal,lengthRes);
% a = 0.0254/16;
% 
% rhoGrid = ones(lengthRes,1)*rhoSpace;
% zGrid = (ones(radialRes,1)*zSpace)';
% 
% M = 2e6;
% C = -mu0/(4*pi)*M;
% 
% 
% for i = 1:lengthRes
%     for j = 1:radialRes
%         
%         rho = rhoSpace(j);
%         z = zSpace(i);
%         
%         integralBz = @(R,Phi) R.*(L/2-z)./(R.^2+(L/2-z)^2+rho^2-2*R.*rho.*cos(Phi)).^(3/2)+...
%             R.*(L/2+z)./(R.^2+(L/2+z)^2+rho^2-2*R.*rho.*cos(Phi)).^(3/2);
%         integralBrho = @(R,Phi) -R.*(2*rho-2*R.*cos(Phi))./(2*(R.^2+(L/2-z).^2+rho^2-2*R*rho.*cos(Phi)).^(3/2))+...
%             R.*(2*rho-2*R.*cos(Phi))./(2*(R.^2+(L/2+z)^2+rho^2-2*R*rho.*cos(Phi)).^(3/2));
%         
%         BzField(i,j) = C*integral2(integralBz,0,a,0,2*pi);
%         BrhoField(i,j) = C*integral2(integralBrho,0,a,0,2*pi);
%     end
% end

radialRes = 400;
lengthRes = 800;
L = 0.0254/2;       % Length of magnet
rhoMaxVal = 0.1;
zMaxVal = 0.1;
rhoMagSpace = linspace(0,rhoMaxVal,radialRes);
zMagSpace = linspace(-zMaxVal,zMaxVal,lengthRes);
a = 0.0254/16;

rhoGrid = ones(lengthRes,1)*rhoMagSpace;
zGrid = (ones(radialRes,1)*zMagSpace)';

BrhoField = importdata('Brho_400by800_L0p0254over2_rho0p0254over16_rhoMaxLMax0p1.csv');
BzField = importdata('Bz_400by800_L0p0254over2_rho0p0254over16_rhoMaxLMax0p1.csv');

%% Magnetic flowmeter setup
numberOfMagnets = 6;        % Designed for 3, additional changes needed for more
magnetSpacing = 60;         % 60 degree spacing for rotating flowmeter
maxOmegaIterations = 10;
rFlowmeter = 0.055;         % Distance from flowmeter center to magnet center
%halfPipeLength = 0.03;
halfPipeLength = 0.03;
dPipe = 0.0391;     % Pipe diameter in meters
%flowSpacing = [-dPipe - 0.0055 - 2*rFlowmeter];        % Distance from magnet tip to liquid metal
flowSpacing = L/2+[0.005 0.01 0.02 0.04 0.08];
dataPoints = length(flowSpacing);

dr = (dPipe/2)/100;                % Pipe grid radial resolution [m]
dTheta = 360/6;                % Pipe grid theta resolution [degrees]
dz = 0.005;

sigma = 3.1e6;              % Galinstan electrical conductivity
Q = LMX_pumpRPM2flowRate(291);            % Flowrate in m^3/s
APipe = pi*(dPipe/2)^2;     % Pipe area in meters^2
uAvg = Q/APipe;             % Average flow speed
angleRes = 60;              % Rotation resolution


% Create magnet rotation angles
for i = 1:numberOfMagnets
    magnetAngles(i,:) = magnetSpacing*(i-1):angleRes:magnetSpacing+magnetSpacing*(i-1)-angleRes;
end

% for i = 1:numberOfMagnets
%     magnetAngles(i,:) = linspace(-magnetSpacing*(i-1),magnetSpacing-magnetSpacing*(i-1),angleRes);
% end

%% Create a mesh around pipe volume

rSpace = 0:dr:dPipe/2;
thetaSpace = 0:dTheta:360-dTheta;
zSpace = -halfPipeLength:dz:halfPipeLength;

rGridSize = length(rSpace);
thetaGridSize = length(thetaSpace);
zGridSize = length(zSpace);
angleCount = size(magnetAngles,2);

xMesh = zeros(thetaGridSize,rGridSize,zGridSize);
yMesh = zeros(thetaGridSize,rGridSize,zGridSize);
zMesh = zeros(thetaGridSize,rGridSize,zGridSize);

% figure(1)
% hold on
% axis equal
% xlabel('xcoors')
% ylabel('ycoors')
% zlabel('zcoors')


for j = 1:rGridSize
    for i = 1:thetaGridSize
        for k = 1:zGridSize
            xMesh(i,j,k) = zSpace(k);
            yMesh(i,j,k) = dPipe/2+rSpace(j)*cosd(thetaSpace(i));
            zMesh(i,j,k) = rSpace(j)*sind(thetaSpace(i));
%             plot3(xMesh(i,j,k),yMesh(i,j,k),zMesh(i,j,k),'bo')
        end
    end
end
% hold off
% error('life')
%% Solve for torque

Torque = zeros(thetaGridSize,rGridSize,zGridSize,angleRes);

aSpace = [20, 100];

aLength = length(aSpace);

totalTorque = 1;
torqueTol = 1e-4;
totalRequiredIterations = aLength*dataPoints*maxOmegaIterations*angleCount*zGridSize*rGridSize*thetaGridSize;
progressIndex = 0;

for aIndex = 1:aLength
    u0 = Q*(aSpace(aIndex)+2)/(2*pi*aSpace(aIndex)*(dPipe/2)^2);
    for dataIndex = 1:dataPoints
        omegaMin = 0;
        omegaMax = 2;
        volumeTotal = 0;
        omegaTolerance = 0.05;
        omegaPrecision = omegaMax-omegaMin;
        totalTorque = 1;
        absTorque = 1;
        yMesh_ = flowSpacing(dataIndex)+rFlowmeter+yMesh;
        iterationCount = 0;
        while iterationCount<=maxOmegaIterations
            omegaGuess = (omegaMax+omegaMin)/2;
            for n = 1:angleCount
                xMag = -rFlowmeter*sind(magnetAngles(:,n));
                yMag = rFlowmeter*cosd(magnetAngles(:,n));
                zMag = 0*ones(size(magnetAngles,1),1);
                counter = 0;
                xGrid_ = ones(numberOfMagnets,1)*reshape(xMesh,1,[])-xMag*ones(1,length(reshape(xMesh,1,[])));
                yGrid_ = ones(numberOfMagnets,1)*reshape(yMesh_,1,[])-yMag*ones(1,length(reshape(yMesh_,1,[])));
                zGrid_ = ones(numberOfMagnets,1)*reshape(zMesh,1,[])-zMag*ones(1,length(reshape(zMesh,1,[])));
                cosRot1Angle = cosd(-magnetAngles(:,n));
                sinRot1Angle = sind(-magnetAngles(:,n));
                cosRot2Angle = cosd(magnetAngles(:,n));
                sinRot2Angle = sind(magnetAngles(:,n));
                
                xRot = cosRot1Angle.*xGrid_ - sinRot1Angle.*yGrid_;
                yRot = sinRot1Angle.*xGrid_ + cosRot1Angle.*yGrid_;
                zRot = zGrid_;
                
                uPipe = -pipeFlowProfileOptimizer(rSpace,dPipe/2,u0,aSpace(aIndex));
                %uPipe = -u0/2;
                uRel_ = ones(numberOfMagnets,1)*reshape(uPipe - (-omegaGuess*yMesh_),1,[]);
                vRel_ = ones(numberOfMagnets,1)*reshape(omegaGuess*xMesh,1,[]);
                wRel_ = ones(numberOfMagnets,1)*reshape(zeros(size(vRel_)),1,[]);
                
                for k = 1:zGridSize
                    for j = 1:rGridSize
                        for i = 1:thetaGridSize
%                             tic
                            counter = counter+1;
                            
                            %% Figure out how to implement new magnetic field...
                            rho = sqrt(xRot(:,counter).^2+zRot(:,counter).^2);
                            Brho_pt = zeros(6,1);
                            Bz_pt = Brho_pt;
                            for magCounter = 1:numberOfMagnets
                                if rho(magCounter) > rhoMaxVal || abs(yRot(magCounter,counter)) > zMaxVal
                                    Brho_pt(magCounter) = 0;
                                    Bz_pt(magCounter) = 0;
                                else
                                    Brho_pt(magCounter) = interp2(rhoGrid,zGrid,BrhoField,rho(magCounter),yRot(magCounter,counter),'spline');
                                    Bz_pt(magCounter) = interp2(rhoGrid,zGrid,BzField,rho(magCounter),yRot(magCounter,counter),'spline');
 
                                end
                            end
                            BxP = Brho_pt.*xRot(:,counter)./sqrt(xRot(:,counter).^2+zRot(:,counter).^2);
                            ByP = Bz_pt;
                            BzP = Brho_pt.*zRot(:,counter)./sqrt(xRot(:,counter).^2+zRot(:,counter).^2);
                            
                            BxP(isnan(BxP)) = 0;
                            ByP(isnan(ByP)) = 0;
                            BzP(isnan(BzP)) = 0;
                            
                            %% Old code
%                             BxP = mu0/(4*pi)*3*pm*(xRot(:,counter).*yRot(:,counter))./(xRot(:,counter).^2+yRot(:,counter).^2+zRot(:,counter).^2).^(5/2);
%                             ByP = mu0/(4*pi)*pm*(2*yRot(:,counter).^2-zRot(:,counter).^2-xRot(:,counter).^2)./(xRot(:,counter).^2+yRot(:,counter).^2+zRot(:,counter).^2).^(5/2);
%                             BzP = mu0/(4*pi)*3*pm*(zRot(:,counter).*yRot(:,counter))./(xRot(:,counter).^2+yRot(:,counter).^2+zRot(:,counter).^2).^(5/2);
                            
                            Bx = cosRot2Angle.*BxP - sinRot2Angle.*ByP;
                            By = sinRot2Angle.*BxP + cosRot2Angle.*ByP;
                            Bz = BzP;
                            
                            if sum(isnan(Bx)) > 0 || sum(isnan(By)) > 0 || sum(isnan(Bz)) > 0
                                error('NaN BOIIIIIIZ')
                            end
                            
                            BMag = sqrt(Bx.^2+By.^2+Bz.^2);
                            
                            uRel = uRel_(:,counter);
                            vRel = vRel_(:,counter);
                            wRel = wRel_(:,counter);
                            
                            jxBForce = -sigma*cross(cross([uRel,vRel,wRel],[Bx,By,Bz]),[Bx,By,Bz]);
                            
                            Torque(i,j,k,n) = (xMag'*jxBForce(:,2) - yMag'*jxBForce(:,1))*dr*(rSpace(j))*dTheta*pi/180*dz; % Volume size may need fixing
                            
                            %                          if uRel(1)<0 && sum(Torque(i,j,k,:))<0
                            %                             error('life')
                            %                          end
%                             totalRequiredIterations = totalRequiredIterations - 1;
%                             stepTime = toc;
%                             disp('Estimated hours until completion:')
%                             disp(stepTime*totalRequiredIterations/3600)
%                             
                                              %error('life')
                        end
                    end
                end
            end
            disp('Flow spacing')
            disp(flowSpacing(dataIndex))
            disp('Iteration')
            disp(iterationCount)
            disp(uRel)
            
            netTorque = sum(Torque,4);
            
            totalTorque = sum(sum(sum(netTorque)));
            
            disp('Omega value')
            disp(omegaGuess)
            disp('Torque')
            disp(totalTorque)
            if totalTorque < 0
                omegaMax = omegaGuess;
            elseif totalTorque > 0
                omegaMin = omegaGuess;
            else

            end
            omegaPrecision = omegaMax-omegaMin;
            absTorque = abs(totalTorque);
            iterationCount = iterationCount+1;
        end
        
        omegaAnswer(dataIndex,aIndex) = omegaGuess;
        aAnswer(dataIndex,aIndex) = aSpace(aIndex);
        flowSpacingAnswer(dataIndex,aIndex) = flowSpacing(dataIndex);
        
    end
    
end


%end