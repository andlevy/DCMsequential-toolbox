function [pE,pC] = spm_L_com_priors(dipfit,pE,pC)

% defaults
%--------------------------------------------------------------------------
try, model    = dipfit.model;    catch, model    = 'BAC'; end
try, type     = dipfit.type;     catch, type     = 'COM'; end
% try, location = dipfit.location; catch, location = 0;     end
try, pC;                         catch, pC       = [];    end





% contributing states (encoded in J)
%==========================================================================
if ischar(model), mod.source = model; model = mod; end
pE.J = {};
pC.J = {};

for i = 1:numel(model)
    
    switch upper(model(i).source)
        
        case{'ERP','SEP'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,9).^0.*-32;               % 9 states
            pC.J{end + 1} = sparse(1,9);  
            
        case{'CMC','TFM'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,8).^0.*-32;           % 8 states
            pC.J{end + 1} = sparse(1,8); 
            
        case{'LFP'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,13).^0.*-32;              % 13 states
            pC.J{end + 1} = sparse(1,13); 
            
        case{'NMM'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,9).^0.*-32;               % 9 states
            pC.J{end + 1} = sparse(1,9);  
            
        case{'NMDA'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,12).^0.*-32;              % 12 states
            pC.J{end + 1} = sparse(1,12); 
            
        case{'CMM'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,12).^0.*-32;              % 12 states
            pC.J{end + 1} = sparse(1,12); 
            
        case{'CMM_NMDA'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,16).^0.*-32;              % 16 states
            pC.J{end + 1} = sparse(1,16);  
            
        case{'MFM'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,36).^0.*-32;              % 36 (9 + 27)
            pC.J{end + 1} = sparse(1,36);        % states
            
        case{'DEM','NFM'}
            %--------------------------------------------------------------
            pE.J{end + 1} = [];                              % null
            pC.J{end + 1} = [];
            
        case{'BGT'}
            %--------------------------------------------------------------
            %assuming data are from STN
            pE.J{end + 1} = sparse(1,10).^0.*-32;               % 10 states
            pC.J{end + 1} = sparse(1,10);
        
        case{'MMC'}
            %--------------------------------------------------------------
            pE.J{end + 1} = sparse(1,8).^0.*-32;   % 8 states
            pC.J{end + 1} = sparse(1,8);

        case {'BCN'}
            pE.J{end + 1} = sparse(1,15).^0.*-32;   % 8 states
            pE.J{end}(15) = 0;
            pC.J{end + 1} = sparse(1,15,1/32,1,15);

        case{'NULL'}
            %--------------------------------------------------------------
            nx   = size(pE.A,1)/m;
            pE.J{end + 1} = ones(1,nx);                      % nx states
            pC.J{end + 1} = ones(1,nx);
        otherwise
            warndlg('Unknown neural model')
            
    end
    
end

