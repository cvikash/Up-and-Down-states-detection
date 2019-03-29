function [K_x,phi_x] = UP_DOWN_DET_compInstPhaseAmpl(x,method)
switch method
    case 'hilbert'
        x_hil = hilbert(x);
        K_x = abs(x_hil);
        phi_x = angle(x_hil); 
    case 'interpol'
        x_hil = hilbert(x);
        K_x = abs(x_hil);
        phi_x_hilbert = angle(x_hil);
        phi_x = zeros(length(x),1);
        zeroCrossings_neg2pos = find(diff(sign(x))==2)+1;
        zeroCrossings_pos2neg = find(diff(sign(x))==-2)+1;
        dx_dt = [0; diff(x)];
        minima = find(diff(sign(dx_dt))==2)+1;
        maxima = find(diff(sign(dx_dt))==-2)+1;
        phi_x(minima) = -pi;
        minima_1=minima-1;
        [~,i2exclude] = intersect(minima_1,maxima);
        minima_1(i2exclude) = [];
        [~,i2exclude2] = intersect(minima_1,zeroCrossings_neg2pos);
        minima_1(i2exclude2) = [];
        [~,i2exclude3] = intersect(minima_1,zeroCrossings_pos2neg);
        minima_1(i2exclude3) = [];
        phi_x(minima_1)=pi;
        phi_x(maxima) = 0;
        phi_x(zeroCrossings_neg2pos)=-pi/2;
        phi_x(zeroCrossings_pos2neg)=pi/2;
        x_specialPoints = sort([minima; minima_1; maxima; zeroCrossings_neg2pos; zeroCrossings_pos2neg]);
        if (x_specialPoints(1)~=1)
            x_specialPoints = [1; x_specialPoints];
            phi_x(1) = phi_x_hilbert(1);
        end
        if (x_specialPoints(end)~=length(x))
            x_specialPoints = [x_specialPoints; length(x)];
            phi_x(end) = phi_x_hilbert(end);
        end        
        phi_x = interp1(x_specialPoints,phi_x(x_specialPoints),[1:1:length(x)])';        
    otherwise
        K_x = [];
        phi_x = [];
end