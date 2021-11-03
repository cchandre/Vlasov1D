function VP1D4f_AnalyzeData
%%
%% Last modified by Cristel Chandre (November 3, 2021)
%% Comments? cristel.chandre@cnrs.fr 
%%
[filename, path] = uigetfile('*.mat');
load([path filename],'suppl','output_modes','Lx','kappa')

P0 = 3;
rho0 = 1;
Nx = length(suppl(1,:))-1;
d = suppl(2,1)-suppl(1,1);
suppl = suppl(1:2*floor(length(suppl(:,1))/2),:);
Nt = length(suppl(:,1));

if strcmp(output_modes,'real')==true
    figure
    imagesc(suppl(:,1),linspace(-Lx,Lx,Nx),suppl(:,2:end).')
    shading flat
    colorbar
    xlabel('$\omega_p t$','interpreter','latex','FontSize',26)
    ylabel('$x$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2,'YDir','normal')    
    figure
    fft_E = fft2(suppl(:,2:end).');
    kx = pi/Lx*(-Nx/2:Nx/2-1);
    om = 2*pi*(-Nt/2:Nt/2-1)/(Nt*d);
    %imagesc(om,kx,log(abs(fftshift(fft_E))))
    imagesc(om,kx,abs(fftshift(fft_E)))
    colorbar
    xlim(2*pi*[-10 10]/(Nt*d))
    ylim([-10 10]*pi/Lx)
    %caxis([-15 log(max(abs(fft_E(:))))])
    xlabel('$\omega$','interpreter','latex','FontSize',26)
    ylabel('$k_x$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2,'YDir','normal') 
    hold on
    lam02 = 1+3*double(kappa)*P0^(2/3)/rho0*kx.^2;
    lambg2 = 1+3*P0/rho0^2*kx.^2;
    lam = 1/sqrt(2).*sqrt(lam02-sqrt(lam02.^2-4*(lam02-lambg2)));
    plot(lam, kx, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    plot(-lam, kx, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    lam = 1/sqrt(2).*sqrt(lam02+sqrt(lam02.^2-4*(lam02-lambg2)));
    plot(lam, kx, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    plot(-lam, kx, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
elseif isinteger(output_modes)==true
    for it=3:output_modes-1
        semilogy(suppl(:,1),abs(suppl(:,it)),'LineWidth',2,'DisplayName',strcat('$k=$ ',' ',num2str(it-2)))
        hold on
    end
    legend('interpreter','latex')
    xlabel('$\omega_p t$','interpreter','latex','FontSize',26)
    ylabel('$E_k(t)$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2)    
end
