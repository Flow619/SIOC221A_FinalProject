function [ACC,VEL,DISP] = compute_DISP_high_pass(ACC,fe,fcut,order)
% function taken from Luke Colosi, Scripps Institution of Oceanography

    Ordercut=order;
    Wcut=2*fcut/fe; 
    [Bcut,Acut] = butter(Ordercut,Wcut,'high');

    Ax=detrend(ACC(1,:));
    Ay=detrend(ACC(2,:));
    Az=detrend(ACC(3,:));

    Axf=detrend(filtfilt(Bcut,Acut,Ax));
    Ayf=detrend(filtfilt(Bcut,Acut,Ay));
    Azf=detrend(filtfilt(Bcut,Acut,Az));

    Vx=cumtrapz(Ax)/fe;
    Vy=cumtrapz(Ay)/fe;
    Vz=cumtrapz(Az)/fe;

    Vxf=detrend(filtfilt(Bcut,Acut,Vx));
    Vyf=detrend(filtfilt(Bcut,Acut,Vy));
    Vzf=detrend(filtfilt(Bcut,Acut,Vz));

    Dx=cumtrapz(Vxf)/fe;
    Dy=cumtrapz(Vyf)/fe;
    Dz=cumtrapz(Vzf)/fe;

    Dxf=detrend(filtfilt(Bcut,Acut,Dx));
    Dyf=detrend(filtfilt(Bcut,Acut,Dy));
    Dzf=detrend(filtfilt(Bcut,Acut,Dz));

    ACC=[Axf;Ayf;Azf];
    VEL=[Vxf;Vyf;Vzf];
    DISP=[Dxf;Dyf;Dzf];

end