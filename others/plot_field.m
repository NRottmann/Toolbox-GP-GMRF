function [] = plot_field(Field,ll,range,max_value,min_value)
% Generates a contour plot of the given field
% Input:
% Field: Field which will be plotted as a counter plot
% ll: number of counter lines
% range: range of the field (D x 1)
% max_value: maximum value for colorbar
% min_value: minimum value for colorbar

[~,XX]=meshgrid(linspace(0,range(1),length(Field(1,:))));
[YY,~]=meshgrid(linspace(0,range(2),length(Field(:,1)))); 
X_Vec = XX(:,1);
Y_Vec = YY(1,:);
contourf(X_Vec,Y_Vec,Field,ll,'LineStyle','none')
caxis([min_value max_value])
colormap(jet)
colorbar

end

