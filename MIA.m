% Variable list:
%
% cell_type_genes: logical matrix in the structure of N rows x M columns, where N = number of
%                  genes in the expression matrix, and M = number of cell
%                  types. There is a '1' in the matrix if that gene (row
%                  index) is a marker for that cell type (column index).
%                  All other matrix elements are '0'. Must be the same
%                  number of genes (N) between scRNA-Seq and ST expression
%                  matrices. 
%
% ST_region_genes: logical matrix in the structure of N rows x P columns, where N =
%                  number of genes in the expression, P = number of tissue
%                  regions. There is a '1' in the matrix if that gene (row
%                  index) is a marker for that tissue region (column
%                  index). All other matrix elements are '0'. 
%
% gene_names: cell vector of length 'N' that contains all gene names from 
%             both scRNA-Seq and ST expression matrices. The script below
%             assumes that the number of genes in both scRNA-Seq and ST are
%             the same. 
%
% cell_type_labels: cell vector of cell type annotations.
%
% tissue_region_labels: cell vector of tissue region annotations.
%
% The output of the code below is a 'MIA' map of enrichment / depletion
% of cell types (rows) in each ST tissue region (columns).

load MIA_example.mat cell_type_genes ST_region_genes gene_names cell_type_labels tissue_region_labels

M = []; M_enr = []; M_dep = [];

cell_type_genes = logical(cell_type_genes);%ST_region_genes对应矩阵，行为基因名，列为spital区域名
ST_region_genes = logical(ST_region_genes);%ST_region_genes对应矩阵，行为基因名，列为spital区域名

figure;
for region = 1 : size(ST_region_genes,2)
    %region = 1
    
    clear ST_gene_list_idx

    ST_gene_list_idx = ST_region_genes(:,region); %取st_gene_list_idx中的第region列(0,1)

    G = gene_names(ST_gene_list_idx); %取ST_gene_list_idx中的为真的基因名
        
    b = length(gene_names); % Number of genes overall / in the gene universe. 转录组总的基因数
    
    for cell_type = 1:size(cell_type_genes,2)
        %cell_type = 1 
        
        C = gene_names(logical(cell_type_genes(:,cell_type)));%取cell_type_genes中的为真的基因名
        a = length(intersect(G,C)); %intersect of spatial and indrop modules单细胞和空间模块基因的交集数
        c = length(G); %number of genes in the spatial module空间模块的基因数
        d = length(C); %number of genes in the combined indrop sets单细胞模块的基因数
        
        M_enr(region,cell_type) = -1*log10(hygecdf(c-a-1,b,c,b-d));    M(region,cell_type) = M_enr(region,cell_type);
        M_dep(region,cell_type) = -1*log10(1-(hygecdf(c-a-1,b,c,b-d)));
        if (M(region,cell_type) < M_dep(region,cell_type))
            M(region,cell_type) = -1*M_dep(region,cell_type);
        end
    end

    
    
    imagesc((M'),[-10 10]); colorbar; %caxis([0 10])
    set(gca,'ytick',1:length(cell_type_labels),'yticklabel',upper(cell_type_labels));
    set(gca,'xtick',1:length(tissue_region_labels),'xticklabel',upper(tissue_region_labels))
    set(gca,'TickLength',[0,0]); set(gca,'fontsize',12)
    yi = [1:length(cell_type_labels)] + 0.5;
    line([0;length(yi)+0.5],[yi;yi], 'color','k','linewidth',1); line([yi;yi], [0;length(yi)+0.5],'color','k','linewidth',1);    
    cb = colorbar; cb.TickDirection = 'out';
        
end