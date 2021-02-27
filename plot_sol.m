function plot_sol(mesh, U)
% AUTEUR : Diallo Amadou, 28/09/2020

max_n_vertices = max(cellfun(@length, mesh.elements));
padding_function = @(E) [E' ...
			NaN(1,max_n_vertices-length(E))];
elements = cellfun(padding_function, mesh.elements, 'UniformOutput', false);
elements = vertcat(elements{:});
Data = [mesh.vertices, U];
patch('Faces', elements, 'Vertices', Data, 'FaceColor', 'interp', 'CData', U );
%axis('square')
xlim([min(mesh.vertices(:,1)) - 0.1, max(mesh.vertices(:,1)) + 0.1])
ylim([min(mesh.vertices(:,2)) - 0.1, max(mesh.vertices(:,2)) + 0.1])
zlim([min(U) - 0.1, max(U) + 0.1])
xlabel('x'); ylabel('y'); zlabel('U');
colormap('jet');
colorbar

end