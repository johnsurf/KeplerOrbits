    function [state_vector, CovSV] = PolyCorrExtrapToFitTime(Iorder, Time, TFit, fit, CovFit, units)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %elements         = zeros(3,rows);
        %derivElements    = zeros(3,rows);
        %secDerivElements = zeros(3,rows);
        elements      = zeros(3,1);
        derivElements = zeros(3,1);
        secDerivElements = zeros(3,1);
        state_vector_sat = [];
        %tIndex = find(track_data(:,2) == tFit);
        %for ii = 1:rows
            ii = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Tinc     = (Time - TFit)/units.TU;
            %Tinc     = track_data(ii,2) - tFit;
            X  = [];
            X1 = [zeros(3)];
            X2 = [zeros(3), zeros(3)];
            Column     = eye(3);
            Column1    = eye(3);
            Column2    = eye(3);
            for j      = 1:Iorder
                X      = [X, Column];
                Column = Column*Tinc;  
            end
            %timeFac  = ones(1,rows);
            %derivFac = ones(1,rows);
            timeFac     = 1;
            derivFac    = 1;
            secDerivFac = 1;
            for i = 1:size(fit,2)
                %elements = elements + fit(:,i)*Tinc.^(i-1);
                deriv1st    = i-1;
                deriv2nd    = (i-2)*deriv1st;
                elements(:,ii) = elements(:,ii) + fit(:,i)*timeFac;
                timeFac        = timeFac.*Tinc;
                if i > 1
                    X1 = [X1, Column1*deriv1st*derivFac];
                    derivElements(:,ii) = derivElements(:,ii) + fit(:,i)*deriv1st*derivFac;
                    derivFac      = derivFac.*Tinc;
                end
                if i > 2
                    X2 = [X2, Column2*deriv2nd*secDerivFac];
                    secDerivElements(:,ii) = secDerivElements(:,ii) + fit(:,i)*deriv2nd*secDerivFac;
                    secDerivFac      = secDerivFac.*Tinc;
                end
            end
        %end
        XSV = [X; X1; X2];
        CovSV = XSV*CovFit*XSV';
        
        %%%%  Satellite 9-State Vector at the median time
        %state_vector_sat = [elements/DU; derivElements/VU; secDerivElements/AU]
        %state_vector = [elements; derivElements/units.TU; secDerivElements/units.TU^2];
        state_vector = [elements; derivElements; secDerivElements];
        state_vector = state_vector./units.X_PVARef';
        state_vector = [state_vector; Time];
        CovSV = CovSV./units.COV9Ref;
    end