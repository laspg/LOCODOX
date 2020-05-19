% DOXY_interp_WOA_2_argo interpolate WOA to argo.
%
% SYNTAX
% [WOA] = DOXY_interp_WOA_2_argo(WOA,argoWork, argo, PorD, nearShore, datat)
%
% DESCRIPTION
% DOXY_interp_WOA_2_argo  interpolate WOA data over argo horizontal and
% vertical grid. The vertical grid could be pressure (interpP) and over
% density (interpD).
%
% INPUT
%
%   WOA (structure)       filled with the variable and its attributes found
%                         in the NetCDF WOA climatology file
%                         Example:
%                         WOA = 
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                 psatwoa: [1x1 struct]
%                                psal_woa: [1x1 struct]
%                                temp_woa: [1x1 struct]
%                                   ..........
%
%   argo (structure)      float data structure
%                         Example:
%                         argo = 
%                                     dimorder: 'C'
%                                    data_type: [1x1 struct]
%                               format_version: [1x1 struct]
%                             handbook_version: [1x1 struct]
%                          reference_date_time: [1x1 struct]
%                                           ...
%                                         juld: [1x1 struct]
%                                      juld_qc: [1x1 struct]
%                                juld_location: [1x1 struct]
%                                     latitude: [1x1 struct]
%                                    longitude: [1x1 struct]
%                                         pres: [1x1 struct]
%                                   molar_doxy: [1x1 struct]
%                                molar_doxy_qc: [1x1 struct]
%                                         doxy: [1x1 struct]
%                                      doxy_qc: [1x1 struct]
%                                doxy_adjusted: [1x1 struct]
%                             doxy_adjusted_qc: [1x1 struct]
%                          doxy_adjusted_error: [1x1 struct]
%                          history_institution: [1x1 struct]
%                                 history_step: [1x1 struct]
%                                            ...
%
%   argowork (structure)  float intermediate working structure
%                         Example:
%                         argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%
%   PorD (cell string)    cell array of string, with the fields name
%                         needed for interpolation
%                         Example: {'preswoa','pres_adjusted'}
%
%   nearShore (boolean)   boolean array indicating if profil is near shore
%                         (1) or not (0)
%
%   datat (array)         numero of the day in the year, for each cycle
%                         numero of a WMO.
%
% OUTPUT
%   WOA (structure)       filled with the variable and its attributes found
%                         in the NetCDF WOA climatology file, and completed
%                         with interpolated WOA data over argo levels
%                         Example:
%                         WOA = 
%                                 dimorder: 'C'
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%                                    depth: [1x1 struct]
%                                     time: [1x1 struct]
%                                  doxywoa: [1x1 struct]
%                                  psatwoa: [1x1 struct]
%                                 psal_woa: [1x1 struct]
%                                 temp_woa: [1x1 struct]
%                                  density: [1x1 struct]
%                                  preswoa: [1x1 struct]
%                                      obj: 'ObsInSitu'
%                                fillisnan: 0
%                              an_dens_woa: [1x1 struct]
%                               doxywoa_CV: [1x1 struct]
%                          psatwoa_interpP: [1x1 struct]
%                       doxywoa_CV_interpP: [1x1 struct]
%                         temp_woa_interpP: [1x1 struct]
%                         psal_woa_interpP: [1x1 struct]
%                          preswoa_interpP: [1x1 struct]
%                        psatwoa_H_interpP: [1x1 struct]
%                     doxywoa_CV_H_interpP: [1x1 struct]
%                       temp_woa_H_interpP: [1x1 struct]
%                       psal_woa_H_interpP: [1x1 struct]
%                        preswoa_H_interpP: [1x1 struct]
%                          psatwoa_interpD: [1x1 struct]
%                       doxywoa_CV_interpD: [1x1 struct]
%                         temp_woa_interpD: [1x1 struct]
%                         psal_woa_interpD: [1x1 struct]
%                          density_interpD: [1x1 struct]
%                        psatwoa_H_interpD: [1x1 struct]
%                     doxywoa_CV_H_interpD: [1x1 struct]
%                       temp_woa_H_interpD: [1x1 struct]
%                       psal_woa_H_interpD: [1x1 struct]
%                        density_H_interpD: [1x1 struct]
%
% CALL
%
% SEE ALSO
% DOXY_corr_main

% HISTORY
%   $created: 03/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:


function [WOA] = DOXY_interp_WOA_2_argo(WOA,argoWork, argo, PorD, nearShore, datat)

% =========================================================================
%% Initialisation
% =========================================================================

if any(~cellfun('isempty',strfind(PorD, 'pres')))
    letterPorD = 'P';
elseif any(~cellfun('isempty',strfind(PorD, 'dens')))
    letterPorD = 'D';
else
     letterPorD = 'P';   
end

VarInWoa = {'psatwoa','doxywoa_CV','temp_woa','psal_woa',PorD{1}};
VarInArgo = {'psat','doxy_adjusted','temp_adjusted','psal_adjusted',PorD{2}};

% ---------------------------------------------------------------------
% Initialize the WOA interpolated variables
% ---------------------------------------------------------------------
cmpl = sprintf('_interp%s',letterPorD);
for v = 1:length(VarInWoa)
    WOA.([VarInWoa{v} cmpl]) = argoWork.(VarInArgo{v});
    WOA.([VarInWoa{v} cmpl]).data = NaN(size(argoWork.(VarInArgo{v}).data));
end

% =====================================================================
%% horizontal interpolation : linear if away from shore, nearest if not
% =====================================================================
nlevels = length(WOA.depth.data);
X = double(WOA.longitude.data);
Y = double(WOA.latitude.data);
T = double(WOA.time.data);

%Extrapolation de T car pèriodique 
T_extr=zeros(length(T)+2,1);
T_extr(1)=T(end)-365.25;
T_extr(end)=T(1)+365.25;
T_extr(2:end-1)=T;

%Extrapolation de la matrice 
for v=1:length(VarInWoa)
    WOA_extr.(VarInWoa{v}).data=zeros([(length(T)+2)  size(squeeze(WOA.(VarInWoa{v}).data(1,:,:,:)))]);
    WOA_extr.(VarInWoa{v}).data(1,:,:,:)=WOA.(VarInWoa{v}).data(end,:,:,:);
    WOA_extr.(VarInWoa{v}).data(end,:,:,:)=WOA.(VarInWoa{v}).data(1,:,:,:);
    WOA_extr.(VarInWoa{v}).data(2:end-1,:,:,:)=WOA.(VarInWoa{v}).data(:,:,:,:);
end

for z=1:nlevels
    
    % distance >1° des cotes    
    Xq = double(argo.longitude.data(logical(~nearShore))) ;
    isneg = Xq < 0;
    if any(isneg)
        Xq(isneg) = Xq(isneg)+360 ;
    end
    Yq=double(argo.latitude.data(logical(~nearShore)));
    Tq=double(datat(logical(~nearShore)));
    
    for v=1:length(VarInWoa)
        WOA.([VarInWoa{v} '_H' cmpl]).data(find(~nearShore),z) = ...
            interpn(T_extr,Y,X,squeeze(WOA_extr.(VarInWoa{v}).data(:,z,:,:)),Tq,Yq,Xq,'linear');
    end
    
    % distance <1° des cotes
    if any(logical(nearShore))
        Xq=double(argo.longitude.data(logical(nearShore)));
        isneg = double(argo.longitude.data(logical(nearShore))) < 0;
        if any(isneg)
            Xq(isneg)=double(argo.longitude.data(logical(nearShore)))+360;
        end
        Yq=double(argo.latitude.data(logical(nearShore)));
        Tq=double(datat(logical(nearShore)));
        
        %ancienne version
        % for v=1:length(VarInWoa)
        %      WOA.([VarInWoa{v} '_H' cmpl]).data(find(nearShore),z) = interpn(T,Y,X,squeeze(WOA.(VarInWoa{v}).data(:,z,:,:)),Tq,Yq,Xq,'nearest');
        % end
        
        %01/07/19 marine
        %Détermination des positions ok, pour lesquelles on a une valeur
        %de pression
        X_ind=[max(find(X<min(Xq))) : min(find(X>max(Xq)))];
        Y_ind=[max(find(Y<min(Yq))) : min(find(Y>max(Yq)))];
        [a,b]=find(isnan(squeeze(WOA.(VarInWoa{5}).data(:,z,Y_ind,X_ind))));
        test_pos=ones(length(Y_ind), length(X_ind));
        test_pos(unique(b))=0;
        pos_ok_X=[];
        pos_ok_Y=[];
        for i=1:length(Y_ind)
            pos_ok_X=[pos_ok_X X_ind(test_pos(i,:)==1)];
            pos_ok_Y=[pos_ok_Y X_ind(test_pos(i,:)==1).*0 + Y_ind(i)];
        end
        
        for v=1:length(VarInWoa)
            if any(any(test_pos))                
                for j=1:length(Tq)
                    vect_dist=dist([Y(pos_ok_Y),Yq(j)],[X(pos_ok_X),Xq(j)]);
                    ind_coord=find(vect_dist==min(vect_dist));

                    
                    WOA2_tmp(j) = interp1(T_extr,WOA_extr.(VarInWoa{v}).data(:,z,pos_ok_Y(ind_coord(1)),pos_ok_X(ind_coord(1))),Tq(j),'linear','extrap');   
%                     end
                    
                    
                    
%                     ind_neg=(T-Tq(j))<0;
%                     diff_time=T-Tq(j);
%                     diff_time(ind_neg)=diff_time(ind_neg)+365.25;
%                     ind_time=find(diff_time==min(diff_time));
%                     WOA2_tmp(j)=WOA.(VarInWoa{v}).data(ind_time,z,pos_ok_Y(ind_coord(1)),pos_ok_X(ind_coord(1)));
                end
                WOA.([VarInWoa{v} '_H' cmpl]).data(logical(nearShore),z)=WOA2_tmp;
            else
                WOA.([VarInWoa{v} '_H' cmpl]).data(logical(nearShore),z)=NaN;
            end
                
        end
              
    end
end

% =====================================================================
%% Vertical interpolation over pressure or density
% =====================================================================
for i=1:size(argoWork.pres_adjusted.data,1)
    isok = ~isnan(WOA.([PorD{1} '_H' cmpl]).data(i,:));
    if any(isok)
        Z = abs(WOA.([PorD{1} '_H' cmpl]).data(i,isok));
    else
        Z = NaN;
    end
    
    if any(isnan(Z))
        WOA.(['psatwoa' cmpl]).data(i,:) = NaN;
        WOA.(['doxywoa_CV' cmpl]).data(i,:) = NaN;
    else
        Zq = argoWork.(PorD{2}).data(i,:);
        for v=1:length(VarInWoa)
            WOA.([VarInWoa{v} cmpl]).data(i,:) = ...
                    interp1(Z',WOA.([VarInWoa{v} '_H' cmpl]).data(i,isok),Zq,'linear');      
        end
    end
end

