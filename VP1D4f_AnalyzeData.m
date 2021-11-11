function VP1D4f_AnalyzeData
%%
%% Last modified by Cristel Chandre (November 4, 2021)
%% Comments? cristel.chandre@cnrs.fr 
%%

[filename, path] = uigetfile('*.mat');
load([path filename],'data','output_modes','Lx','kappa','precision')

P0 = 1;
rho0 = 1;
Nx = length(data(1,:))-1;
d = data(2,1)-data(1,1);
data = data(1:2*floor(length(data(:,1))/2),:);
Nt = length(data(:,1));

if strcmp(output_modes,'real')==true
    figure
    imagesc(data(:,1),linspace(-Lx,Lx,Nx),data(:,2:end).')
    shading flat
    cmap1 = 1-hot(256);
    cmap2 = flipud(hot(256));
    colormap([flip(cmap1(2:end,:));cmap2])
    colorbar
    xlabel('$\omega_p t$','interpreter','latex','FontSize',26)
    ylabel('$x$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2,'YDir','normal')    
    figure
    fft_data = fftshift(fft2(data(:,2:end).'));
    fft_data(abs(fft_data)<=precision*max(abs(fft_data(:)))) = NaN;
    kx = pi/Lx*(-Nx/2:Nx/2-1);
    om = 2*pi*(-Nt/2:Nt/2-1)/(Nt*d);
    imagesc(om,kx,log10(abs(fft_data)))
    colormap(flipud(pink))
    colorbar
    xlim(2*pi*[-200 200]/(Nt*d))
    ylim([-10 10]*pi/Lx)
    caxis([log10(precision) log(max(abs(fft_data(:))))])
    xlabel('$\omega$','interpreter','latex','FontSize',26)
    ylabel('$k_x$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2,'YDir','normal') 
    hold on
    lam02 = 1+3*double(kappa)*P0^(2/3)/rho0*kx.^2;
    lambg2 = 1+3*P0/rho0^2*kx.^2;
    lam = sqrt(lam02-sqrt(lam02.^2-4*(lam02-lambg2)))/sqrt(2);
    plot(lam(imag(lam)==0), kx(imag(lam)==0), 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    plot(-lam(imag(lam)==0), kx(imag(lam)==0), 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    lam = sqrt(lam02+sqrt(lam02.^2-4*(lam02-lambg2)))/sqrt(2);
    plot(lam(imag(lam)==0), kx(imag(lam)==0), 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
    plot(-lam(imag(lam)==0), kx(imag(lam)==0), 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
elseif isinteger(output_modes)==true
    for it=3:output_modes-1
        semilogy(data(:,1),abs(data(:,it)),'LineWidth',2,'DisplayName',strcat('$k=$ ',' ',num2str(it-2)))
        hold on
    end
    legend('interpreter','latex')
    xlabel('$\omega_p t$','interpreter','latex','FontSize',26)
    set(gca,'box','on','FontSize',20,'LineWidth',2)    
end
