function [J] = nearby_correct_kop(m,I_non_correct)
% function [J] = Nearby_Correct_Kop(m,I_non_correct,I_ST,I_EAc,I_EAo,I_VAc,I_VAo)
%%
J = zeros(6,6,3);
%
ST  = 1;
EAc = 2;
EAo = 3;
VAc = 4;
VAo = 5;
%
% ind_near = 0;
if m == [1,1,1]; 
    
%     ind_near = 1;
    
    C = [ ST  0.0 EAo EAo EAo EAo
          0.0 ST  EAo EAo EAo EAo
          EAo EAo ST  0.0 EAo EAo
          EAo EAo 0.0 ST  EAo EAo
          EAo EAo EAo EAo ST  0.0
          EAo EAo EAo EAo 0.0 ST  ];
      
elseif m == [2,1,1]; 
    
%     ind_near = 1;
    
    C = [ 0.0 ST  EAo EAo EAo EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 EAo EAc 0.0 VAo VAo
          0.0 EAo 0.0 EAc VAo VAo
          0.0 EAo VAo VAo EAc 0.0
          0.0 EAo VAo VAo 0.0 EAc ];
      
 elseif m == [2,1,2]; 
    
%      ind_near = 1;
     
    C = [ 0.0 EAc VAo VAo 0.0 EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo VAc 0.0 0.0 VAo
          0.0 VAo 0.0 VAc 0.0 VAo
          0.0 EAo VAo VAo 0.0 EAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];
      
 elseif m == [1,1,2]; 
    
%      ind_near = 1;
     
    C = [ EAc 0.0 VAo VAo 0.0 EAo
          0.0 EAc VAo VAo 0.0 EAo
          VAo VAo EAc 0.0 0.0 EAo
          VAo VAo 0.0 EAc 0.0 EAo
          EAo EAo EAo EAo 0.0 ST
          0.0 0.0 0.0 0.0 0.0 0.0 ];
      
 elseif m == [1,2,1]; 
    
%      ind_near = 1;
     
    C = [ EAc 0.0 0.0 EAo VAo VAo
          0.0 EAc 0.0 EAo VAo VAo
          EAo EAo 0.0 ST  EAo EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          VAo VAo 0.0 EAo EAc 0.0
          VAo VAo 0.0 EAo 0.0 EAc ];

elseif m == [2,2,1]; 
    
%      ind_near = 1;
     
    C = [ 0.0 EAc 0.0 EAo VAo VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 EAo 0.0 EAc VAo VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAo VAc 0.0
          0.0 VAo 0.0 VAo 0.0 VAc ];

 elseif m == [1,2,2]; 
    
%      ind_near = 1;
     
    C = [ VAc 0.0 0.0 VAo 0.0 VAo
          0.0 VAc 0.0 VAo 0.0 VAo
          VAo VAo 0.0 EAc 0.0 EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          VAo VAo 0.0 EAo 0.0 EAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];

 elseif m == [2,2,2]; 
    
%      ind_near = 1;
     
    C = [ 0.0 VAc 0.0 VAo 0.0 VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAc 0.0 VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAo 0.0 VAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ind_near == 1
    for ii = 1 : 6
        for jj = 1 : 6
            ind = C(ii,jj);
            %         switch ind
            %             case 0
            %                 J(ii,jj,:) = I_non_correct(ii,jj,:);
            %             case 1
            %                 J(ii,jj,:) = I_ST(:,1);
            %             case 2
            %                 J(ii,jj,:) = I_EAc(:,1);
            %             case 3
            %                 J(ii,jj,:) = I_EAo(:,1);
            %             case 4
            %                 J(ii,jj,:) = I_VAc(:,1);
            %             case 5
            %                 J(ii,jj,:) = I_VAo(:,1);
            %         end
            switch ind
                case 1
                    J(ii,jj,:) = 0.0;
                otherwise
                    J(ii,jj,:) = I_non_correct(ii,jj,:);
            end
            
        end
    end
% else
%     J = I_non_correct;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%