function [phi pt] = genphi(N1,N2, c)
% c: Number of coils
phi = zeros(N1, N2, c);
phi0 = -pi + 2*pi*rand(c, 6);     %1 0.5 to 2 1
phi0(:,2:3) = phi0(:,2:3).*1.8/N1;
phi0(:,4:6) = phi0(:,4:6).*0.9/N1/N2;

% phi0 = -3*pi + 6*pi*rand(c, 6);
% phi0(:,2:3) = phi0(:,2:3).*1/N;
% phi0(:,4:6) = phi0(:,4:6).*0.5/N/N;

%  phi0(:,1) = 0;
for i = 1:N1,
		for j = 1:N2,
            for k = 1:c,
                  phi(i, j, k) = phi0(k,1) + phi0(k, 2)*(i-N1/2) + phi0(k, 3)*(j-N2/2)...
                                +phi0(k, 4)*((i-N1/2)^2) + phi0(k, 5)*((j-N2/2)^2)...
                                +phi0(k, 6)*((i-N1/2)*(j-N2/2));
%                     phi(i, j, k) = 0;
            end
        end
end
pt = phi0(:,1);
end