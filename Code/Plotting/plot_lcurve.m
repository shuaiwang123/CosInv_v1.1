function lambda = plot_lcurve(missfit,roughness,lap,UNIT2cm)
% Using the L-curve method to dertimine the Best Laplacian
% 
% MISSFIT   sum of square of the residuals between model data and observed data
% ROUGHNESS  the roughness of the slip distribution 
% LAP  smoothing factor to balance RMS and ROUGNNESS
%
% ###
%   Note that the original unit of the input variables for MISSFIT and
%   ROUGHNESS are m^2 and m/km^2 respectively. For the beautiful of figure, 
%   we transform the unit of these two variables to cm^2 and cm/km^2
%   repectively.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_lap     = numel(lap);
missfit   = missfit*UNIT2cm^2;
roughness = roughness*UNIT2cm;

figure;
box on;axis on;
plot(roughness,missfit,'-bs','MarkerFaceColor','k','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',8);
for p = 1:2:n_lap;
    text(roughness(p),missfit(p),['  k =',num2str(lap(p))],'color','b');
end
xlabel('Roughness(cm/km^2)');ylabel('Missfit(cm^2)');title('Choose the best Laplacian weight through L-curve');
hold on;

id = input('    [Choose the ID number of the best Laplacian weight]');
lambda = lap(id);
plot(roughness(id),missfit(id),'-bd','MarkerFaceColor','k','LineWidth',2,'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',10);
text(roughness(id),missfit(id),['  k =',num2str(lap(id))],'color','b');
hold off;
disp(['      ... The best Laplacian weight is ' num2str(lambda) ' ...']);