function RLFF_InversionFunction_041018
clear
clc
close all

aGuess = 5;
omegaData = [1.2 1.3 1.4 1.5 1.6];

pm = 10;                    % Magic magnet number
mu0 = 4*pi*1e-7;            % Permittivity of free space
BrhoField = importdata('ActualMagnetDims_Brho.csv');
BzField = importdata('ActualMagnetDims_Bz.csv');
radialRes = 400;
lengthRes = 800;
L = 0.0254/2;       % Length of magnet
rhoMaxVal = 0.1;
zMaxVal = 0.1;
rhoSpace = linspace(0,rhoMaxVal,radialRes);
zSpace = linspace(-zMaxVal,zMaxVal,lengthRes);

flowSpacing = [L/2+0.00];

rhoGrid = ones(lengthRes,1)*rhoSpace;
zGrid = (ones(radialRes,1)*zSpace)';

a = 0.0254*3/16/2;
numberOfMagnets = 6;        % Designed for 3, additional changes needed for more
magnetSpacing = 60;         % 60 degree spacing for rotating flowmeter
rFlowmeter = 0.055;         % Distance from flowmeter center to magnet center
halfPipeLength = 0;
dPipe = 0.0391;     % Pipe diameter in meters
dataPoints = length(flowSpacing);

dr = (dPipe/2)/64;                % Pipe grid radial resolution [m]
dTheta = 360/2;                % Pipe grid theta resolution [degrees]
dz = 0.03;

sigma = 3.1e6;              % Galinstan electrical conductivity
Q = LMX_pumpRPM2flowRate(291);            % Flowrate in m^3/s
APipe = pi*(dPipe/2)^2;     % Pipe area in meters^2
wPipe = 0.006;              % Pipe wall thickness
uAvg = Q/APipe;             % Average flow speed
angleRes = 60;              % Rotation resolution


% Create magnet rotation angles
for i = 1:numberOfMagnets
    magnetAngles(i,:) = magnetSpacing*(i-1):angleRes:magnetSpacing+magnetSpacing*(i-1);
end

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

%% Solve for torque

Torque = zeros(thetaGridSize,rGridSize,zGridSize,angleRes);

omegaLength = length(omegaData);

totalTorque = 1;
maxOmegaIterations = 10;
totalRequiredIterations = omegaLength*dataPoints*maxOmegaIterations*angleCount*zGridSize*rGridSize*thetaGridSize;
progressIndex = 0;

for omegaIter = 1:omegaLength
aActual(omegaIter) = fminsearch(@solveTorque,aGuess);

disp(aActual)
end





    function totalTorqueCost = solveTorque(aGuess)
        
        u0 = Q*(aGuess+2)/(2*pi*aGuess*(dPipe/2)^2);
        totalTorque = 1;
        yMesh_ = flowSpacing+wPipe+rFlowmeter+yMesh;
        omegaGuess = omegaData(omegaIter);
        for n = 1:angleCount
            xMag = -(rFlowmeter-L/2)*sind(magnetAngles(:,n));
            yMag = (rFlowmeter-L/2)*cosd(magnetAngles(:,n));
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
            
            uPipe = -pipeFlowProfileOptimizer(rSpace,dPipe/2,u0,aGuess);
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
                        
                        Bx = cosRot2Angle.*BxP - sinRot2Angle.*ByP;
                        By = sinRot2Angle.*BxP + cosRot2Angle.*ByP;
                        Bz = BzP;
                        
                        if sum(isnan(Bx)) > 0 || sum(isnan(By)) > 0 || sum(isnan(Bz)) > 0
                            error('NaN')
                        end
                        
                        uRel = uRel_(:,counter);
                        vRel = vRel_(:,counter);
                        wRel = wRel_(:,counter);
                        
                        jxBForce = -sigma*cross(cross([uRel,vRel,wRel],[Bx,By,Bz]),[Bx,By,Bz]);
                        
                        Torque(i,j,k,n) = ( (xMag'*jxBForce(:,2) - yMag'*jxBForce(:,1)) )*dr*(rSpace(j))*dTheta*pi/180*dz; % Volume size may need fixing
                        
                        %                          if uRel(1)<0 && sum(Torque(i,j,k,:))<0
                        %                             error('life')
                        %                          end
                        %                             totalRequiredIterations = totalRequiredIterations - 1;
                        %                             stepTime = toc;
                        %                             disp('Estimated hours until completion:')
                        %                             disp(stepTime*totalRequiredIterations/3600)
                        %
                    end
                end
            end
            %             disp('Flow spacing')
            %             disp(flowSpacing-L/2)
            %             disp('Iteration')
            %             disp(iterationCount)
        end
        netTorque = sum(Torque,4);
        airDragTorque = 8*pi*rFlowmeter^3*1.81e-5*omegaGuess;
        
        totalTorque = sum(sum(sum(netTorque)))-12*airDragTorque;
        
        totalTorqueCost = totalTorque^2;
        
        fprintf('current a value: %d\n',aGuess)
        fprintf('cost: %d\n\n',totalTorqueCost)
        
        
    end

end