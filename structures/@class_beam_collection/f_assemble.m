%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj=f_assemble(obj,loadstep,eval_nonlin)
            add_eigenmass=obj.settings.gravity;
            add_fuelmass=obj.settings.fuel_mass;
            add_engineforces=obj.settings.engine;
            add_gearforces=obj.settings.landing_gear;
            
            sys_ndof=0;
            updateTotalK=0;
            updateTotalM=0;
            updateTotalQ=0;
            %assemble and set boundary conditions in subsystems first
            for i=1:1:length(obj.beam)
                if eval_nonlin==1
                    obj.beam(i)=obj.beam(i).f_nonlinassemble(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces,loadstep);
                else
                    updateTotalK=(updateTotalK || obj.beam(i).update_K);
                    updateTotalM=(updateTotalM || obj.beam(i).update_M);
                    updateTotalQ=(updateTotalQ || obj.beam(i).update_Q);
                    obj.beam(i)=obj.beam(i).f_assemble(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces);
                end
                if updateTotalK==0 && updateTotalM==0 && updateTotalQ==1
                    obj.beam(i)=obj.beam(i).f_set_BC('onlyQ');
                else
                    obj.beam(i)=obj.beam(i).f_set_BC();
                end
            end
            
            % calculate overall system DOFs
            for i=1:1:length(obj.beam)
                sys_ndof(i)=length(obj.beam(i).Kff);
            end
            
            % preallocate system stiffness and load and mass matrices
            obj.Kglob=zeros(sum(sys_ndof));
            
            obj.Mglob=zeros(sum(sys_ndof));
            obj.Mglob_lumped=zeros(sum(sys_ndof));
            
            obj.Fglob=zeros(sum(sys_ndof),1);
            
            % for nonlinear iteration internal force vector is required as
            % well
            if eval_nonlin==1
                obj.Pglob=zeros(sum(sys_ndof),1);
            end
            
            %
            obj.sort_vec=zeros(sum(sys_ndof),2);
            obj.dof_node_beam=zeros(sum(sys_ndof),3);
            
            sys_pos=cumsum(sys_ndof);
            sys_pos=[0 sys_pos(1:end-1)];
            
            for i=1:1:length(obj.beam)
                obj.Kglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))=obj.Kglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))+obj.beam(i).Kff;
                obj.Mglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))=obj.Mglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))+obj.beam(i).Mff; % for mass matrix
                obj.Mglob_lumped(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))=obj.Mglob_lumped(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),sys_pos(i)+1:sys_pos(i)+sys_ndof(i))+obj.beam(i).Mff_lumped; % for mass matrix
                obj.Fglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)=obj.Fglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)+obj.beam(i).Ff;
                if eval_nonlin==1
                    obj.Pglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)=obj.Pglob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)+obj.beam(i).Pf;
                end
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),3)=i;
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),2)=obj.beam(i).Kff_node_dof_info(:,1);
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)=obj.beam(i).Kff_node_dof_info(:,2);
            end
            
            el_ndof=6;
            
            ctr=1;
            for i=1:1:length(obj.coupling_condition)
                for j=1:1:el_ndof
                    if(obj.coupling_condition(i).dof(j)==1)
                        for k=1:1:length(obj.coupling_condition(i).beam_idx)
                            idx=find(obj.beam(obj.coupling_condition(i).beam_idx(k)).Kff_node_dof_info(:,1)==obj.coupling_condition(i).node_idx(k),1,'first');
                            sort_idx=sys_pos(obj.coupling_condition(i).beam_idx(k))+idx+j-1;
                            obj.sort_vec(sort_idx,1)=ctr;
                            obj.sort_vec(sort_idx,2)=sort_idx;
                        end
                        ctr=ctr+1;
                    end
                end
            end
            
%            ctr=1;
%             for i=1:1:length(obj.coupling_condition)
%                 for j=1:1:el_ndof
%                     for k=1:1:length(obj.coupling_condition(i).beam_idx)
%                         if k<3
%                             dof_idx=1;
%                         else
%                             dof_idx=k-1;
%                         end
%                         if(obj.coupling_condition(i).dof(dof_idx,j)==1)
%                             idx=find(obj.beam(obj.coupling_condition(i).beam_idx(k)).Kff_node_dof_info(:,1)==obj.coupling_condition(i).node_idx(k),1,'first');
%                             sort_idx=sys_pos(obj.coupling_condition(i).beam_idx(k))+idx+j-1;
%                             obj.sort_vec(sort_idx,1)=ctr;
%                             obj.sort_vec(sort_idx,2)=sort_idx;
%                             ctr=ctr+1;
%                         end
%                     end
%                 end
%             end
            
            j=ctr;
            for i=1:sum(sys_ndof)
                if obj.sort_vec(i,1)==0
                    obj.sort_vec(i,1)=j;
                    obj.sort_vec(i,2)=i;
                    j=j+1;
                end
            end
            
            len=sum(sys_ndof);
            col_ndof=max(obj.sort_vec(:,1));
            tot_size=sum(sys_ndof);
            obj.n_dof_sys=tot_size;
            
            if updateTotalK
                Ksort=[obj.Kglob obj.sort_vec(:,1)];
                Ksort=sortrows(Ksort,len+1);
                %!!!
                K=zeros(col_ndof,tot_size+1);
                for i=1:1:col_ndof
                    idx=find(Ksort(:,end)==i);
                    if length(idx)>1
                        K(i,:)=sum(Ksort(idx,1:end));
                        K(i,end)=i;
                    else
                        K(i,:)=Ksort(idx,1:end);
                    end
                end

                %!!!
                Ksort=K(:,1:end-1)';
                Ksort=[Ksort obj.sort_vec(:,1)];

                Ksort=sortrows(Ksort,col_ndof+1);

                K=zeros(col_ndof);
                for i=1:1:col_ndof
                    idx=find(Ksort(:,end)==i);
                    if length(idx)>1
                        K(i,:)=sum(Ksort(idx,1:end-1));
                    else
                        K(i,:)=Ksort(idx,1:end-1);
                    end
                end
                obj.Kff=K';
                obj.Kglob=Ksort;
            end
            
            if updateTotalM
                % for mass matrix
                %mass matrix --> Msort
                Msort=[obj.Mglob obj.sort_vec(:,1)];
                Msort=sortrows(Msort,len+1);

                M=zeros(col_ndof,tot_size+1);
                for i=1:1:col_ndof
                    idx=find(Msort(:,end)==i);
                    if length(idx)>1
                        M(i,:)=sum(Msort(idx,1:end));
                        M(i,end)=i;
                    else
                        M(i,:)=Msort(idx,1:end);
                    end
                end

                Msort=M(:,1:end-1)';
                Msort=[Msort obj.sort_vec(:,1)];

                Msort=sortrows(Msort,col_ndof+1);

                M=zeros(col_ndof);
                for i=1:1:col_ndof
                    idx=find(Msort(:,end)==i);
                    if length(idx)>1
                        M(i,:)=sum(Msort(idx,1:end-1));
                    else
                        M(i,:)=Msort(idx,1:end-1);
                    end
                end


                % for lumped mass matrix
                %mass matrix --> Msort
                Msort_lumped=[obj.Mglob_lumped obj.sort_vec(:,1)];
                Msort=sortrows(Msort_lumped,len+1);

                Mlumped=zeros(col_ndof,tot_size+1);
                for i=1:1:col_ndof
                    idx=find(Msort_lumped(:,end)==i);
                    if length(idx)>1
                        Mlumped(i,:)=sum(Msort_lumped(idx,1:end));
                        Mlumped(i,end)=i;
                    else
                        Mlumped(i,:)=Msort_lumped(idx,1:end);
                    end
                end

                Msort_lumped=Mlumped(:,1:end-1)';
                Msort_lumped=[Msort_lumped obj.sort_vec(:,1)];

                Msort_lumped=sortrows(Msort_lumped,col_ndof+1);

                Mlumped=zeros(col_ndof);
                for i=1:1:col_ndof
                    idx=find(Msort_lumped(:,end)==i);
                    if length(idx)>1
                        Mlumped(i,:)=sum(Msort_lumped(idx,1:end-1));
                    else
                        Mlumped(i,:)=Msort_lumped(idx,1:end-1);
                    end
                end

                obj.Mff=M';
                obj.Mff_lumped= Mlumped;
                obj.Mglob=Msort; 
            end
            if updateTotalQ
                %forces
                Fsort=[obj.Fglob obj.sort_vec(:,1)];
                Fsort=sortrows(Fsort,2);
                F=zeros(col_ndof,1);
                for i=1:1:max(obj.sort_vec(:,1))
                    idx=find(Fsort(:,end)==i);
                    if length(idx)>1
                        F(i,:)=sum(Fsort(idx,1:end-1));
                    else
                        F(i,:)=Fsort(idx,1:end-1);
                    end
                end
                %nonlinear
                if eval_nonlin==1
                    Psort=[obj.Pglob obj.sort_vec(:,1)];
                    Psort=sortrows(Psort,2);
                    P=zeros(col_ndof,1);
                    for i=1:1:max(obj.sort_vec(:,1))
                        idx=find(Psort(:,end)==i);
                        if length(idx)>1
                            P(i,:)=sum(Psort(idx,1:end-1));
                        else
                            P(i,:)=Psort(idx,1:end-1);
                        end
                    end
                    obj.Ptest=P;
                end
                
                obj.Ftest=F;
            end
            
            
            %                 obj.Ktest=K';
            %                 obj.Mtest=M'; % Mass Matrix test
            
            
        end
