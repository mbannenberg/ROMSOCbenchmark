% ------------------------------------------------------------------------------
% Netlist parser function. Output given in quasi linear implicit DAE form.
% 
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [E, A, func_p, func_r, x0] = parse_netlist(fname)
    fileID = fopen(fname);
    tline = fgetl(fileID);
   
    %Initialize
    numR=0;
    numL=0;
    numC=0;
    numV=0;     %Number of independent voltage sources
    numD=0;     %Number of diodes
    numM=0;     %Number of MOSFETs
    numI=0;     %Number of independent current sources
    numT=0;     %Number of transistors
    numNode=1;  %Number of nodes, not including ground (node 0).
    
    
    Nodes = containers.Map(['SUB'],[0]);
    while ischar(tline)
        args = strsplit(tline,' ');
        switch(tline(1))
            case {'R'}
                numR=numR+1;
                Resistor(numR).Name = char(args(1));

                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);

                Resistor(numR).Node1= Nodes(char(args(2)));
                Resistor(numR).Node2= Nodes(char(args(3)));
                Resistor(numR).Value= str2double(args(4));
            case {'L'}
                numL=numL+1;
                Inductor(numL).Name = char(args(1));

                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);

                Inductor(numL).Node1= Nodes(char(args(2)));
                Inductor(numL).Node2= Nodes(char(args(3)));
                Inductor(numL).Value= str2double(args(4));
            case {'C'}
                numC=numC+1;
                Capacitor(numC).Name = char(args(1));

                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);

                Capacitor(numC).Node1= Nodes(char(args(2)));
                Capacitor(numC).Node2= Nodes(char(args(3)));
                Capacitor(numC).Value= str2double(args(4));

            case 'V'
                numV=numV+1;
                Vsource(numV).Name = char(args(1));

                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);

                Vsource(numV).Node1= Nodes(char(args(2)));
                Vsource(numV).Node2= Nodes(char(args(3)));
                Vsource(numV).Type = char(args(4));
                Vsource(numV).Value= str2double(args(5));
                if strcmp(Vsource(numV).Type,'sin')
                    Vsource(numV).Freq = str2double(args(6));
                end

            case 'D'    
                numD=numD+1;
                Diode(numD).Name = char(args(1));

                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);

                Diode(numD).Node1= Nodes(char(args(2)));
                Diode(numD).Node2= Nodes(char(args(3)));
                Diode(numD).Value= str2double(args(4));
            case 'M'
                numM=numM+1;
                Mosfet(numM).Name = char(args(1));


                [Nodes, numNode] = addNode(Nodes,char(args(2)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(3)),numNode);
                [Nodes, numNode] = addNode(Nodes,char(args(4)),numNode);

                % Gate Drain Source
                Mosfet(numM).Node1= Nodes(char(args(2)));
                Mosfet(numM).Node2= Nodes(char(args(3)));
                Mosfet(numM).Node3= Nodes(char(args(4)));
                Mosfet(numM).Type = char(args(5));
                Mosfet(numM).Value= str2double(args(6));
            case 'I'
                numI=numI+1;
                Isource(numI).Name=Name{i};
                Isource(numI).Node1=Nodes(N1{i});
                Isource(numI).Node2=Nodes(N2{i});
                try
                    Isource(numI).Value=str2num(arg3{i});
                catch
                    Isource(numI).Value=nan;
                end
            case 'T'
                numT=numT+1;
        end
        tline = fgetl(fileID);
        if numNode == 199
            a = 3;
        end
    end
    
    % Determine size of the system
    numNode = numNode - 1;
    sizeX = numNode + numL + numV;

    % Allocate initial conditions
    x0 = zeros(sizeX,1);


    %% Construct E matrix
    E = zeros(sizeX);
    A_C = zeros(numNode,numC);
    C_diag = zeros(numC);
    for i=1:numC
        n1=Capacitor(i).Node1;
        n2=Capacitor(i).Node2;
        C_diag(i,i) = Capacitor(i).Value();
        if (n1~=0) && (n2~=0)
            A_C(n1,i)=A_C(n1,i) + 1;
            A_C(n2,i)=A_C(n2,i) - 1;
        elseif (n1~=0)
            A_C(n1,i)= A_C(n1,i) - 1;
        elseif (n2~=0)
            A_C(n2,i)= A_C(n2,i) + 1;
        end 
    end

    A_L = zeros(numNode,numL);
    L_diag = zeros(numL);
    for i=1:numL
        n1=Inductor(i).Node1;
        n2=Inductor(i).Node2;
        L_diag(i,i) = Inductor(i).Value();
        if (n1~=0) && (n2~=0)
            A_L(n1,i)=A_L(n1,i) + 1;
            A_L(n2,i)=A_L(n2,i) - 1;
        elseif (n1~=0)
            A_L(n1,i)= A_L(n1,i) - 1; 
        elseif (n2~=0)
            A_L(n2,i)= A_L(n2,i) + 1;
        end 
    end

    E(1:numNode,1:numNode) = A_C*(C_diag*A_C');
    E(numNode+1:numNode+numL,numNode+1:numNode+numL) = L_diag;

    %% Construct A matrix
    A = zeros(sizeX);

    A_R = zeros(numNode,numR);
    R_diag = zeros(numR);
    for i=1:numR
        n1=Resistor(i).Node1;
        n2=Resistor(i).Node2;
        R_diag(i,i) = 1/Resistor(i).Value();
        if (n1~=0) && (n2~=0)
            A_R(n1,i)=A_R(n1,i) + 1;
            A_R(n2,i)=A_R(n2,i) - 1;
        elseif (n1~=0)
            A_R(n1,i)= A_R(n1,i) - 1;
        elseif(n2~=0)
            A_R(n2,i)= A_R(n2,i) + 1;
        end 
    end

    A_V = zeros(numNode,numV);
    r = cell(sizeX,1);
    for i = 1:sizeX
        r{i} = '@(t) 0';
    end
    for i=1:numV
        n1=Vsource(i).Node1;
        n2=Vsource(i).Node2;

        if (n1~=0) && (n2~=0)
            A_V(n1,i)=A_V(n1,i) + 1;
            A_V(n2,i)=A_V(n2,i) - 1;
        elseif (n1~=0)
            A_V(n1,i)= A_V(n1,i) - 1;
        elseif(n2~=0)
            A_V(n2,i)= A_V(n2,i) + 1;
        end 
        if strcmp(Vsource(i).Type, 'cons')
%             r{numNode+numL+i} = @(t) Vsource(i).Value;
            str = sprintf(' + %e',Vsource(i).Value);
            r{numNode+numL+i,1} = strcat(r{numNode+numL+i,1} , str);
            
            x0(n2) = Vsource(i).Value;
        elseif strcmp(Vsource(i).Type, 'v_in')
            str = sprintf(' + v_in(t)');
            r{numNode+numL+i,1} = strcat(r{numNode+numL+i,1} , str);
            
        elseif strcmp(Vsource(i).Type, 'sin')
            str = sprintf(' + %e*sin(%e*2*pi*t)',Vsource(i).Value,Vsource(i).Freq);
            r{numNode+numL+i,1} = strcat(r{numNode+numL+i,1} , str);
        else
            disp('Wrong type of voltage source');
        end
    end
    
    for i = 1:sizeX
        r{i} = str2func(r{i});
    end
    func_r = r;
%     B = strcat(regexprep(cellfun(@func2str, r, 'uni', 0), '^@\(t\)', ''), ';');
%     func_r = str2func(strcat('@(t) [', B{:}, ']'));

    
    
    A(1:numNode,1:numNode) = A_R*R_diag*A_R';
    A(1:numNode,numNode+1:numNode+numL) = A_L;
    A(1:numNode,numNode+numL+1:numNode+numL+numV) = A_V;

    A(numNode+1:numNode+numL,1:numNode)= -A_L';
    A(numNode+numL+1:numNode+numL+numV,1:numNode)=-A_V';

    p = cell(sizeX,1);
    for i = 1:sizeX
        p{i} = '@(x) 0';
    end

    for i = 1:numD
        n1=Diode(i).Node1;
        n2=Diode(i).Node2;

        if (n1~=0) && (n2~=0)
            str = sprintf(' + diode(%e,x(%d),x(%d))',Diode(i).Value,n1,n2);
            p{n1,1} = strcat(p{n1,1} , str);
            str = sprintf(' - diode(%e,x(%d),x(%d))',Diode(i).Value,n1,n2);
            p{n2,1} = strcat(p{n2,1} , str);
        end

        if (n1==0)
            % Might be wrong
            str = sprintf(' - diode(%e,0,x(%d))',Diode(i).Value,n2);
            p{n2,1} = strcat(p{n2,1} , str);
        end
        if (n2==0)
            str = sprintf(' - diode(%e,x(%d),0)',Diode(i).Value,n1);
            p{n1,1} = strcat(p{n1,1} , str);
        end 

    end

    for i = 1:numM
        n1=Mosfet(i).Node1;
        n2=Mosfet(i).Node2;
        n3=Mosfet(i).Node3;

        if strcmp(Mosfet(i).Type,'N')
            if (n2~=0) && (n3~=0)
                str = sprintf(' + mosfet_n(x(%d),x(%d),x(%d),%e)',n1,n2,n3,Mosfet(i).Value);
                p{n2,1} = strcat(p{n2,1} , str);
                str = sprintf(' - mosfet_n(x(%d),x(%d),x(%d),%e)',n1,n2,n3,Mosfet(i).Value);
                p{n3,1} = strcat(p{n3,1} , str);
            end

            if (n2==0)
                str = sprintf(' - mosfet_n(x(%d),0,x(%d),%e)',n1,n3,Mosfet(i).Value);
                p{n3,1} = strcat(p{n3,1} , str);
            end
            if (n3==0)
                str = sprintf(' + mosfet_n(x(%d),x(%d),0,%e)',n1,n2,Mosfet(i).Value);
                p{n2,1} = strcat(p{n2,1} , str);
            end 
        else
            if (n2~=0) && (n3~=0)
                str = sprintf(' - mosfet_p(x(%d),x(%d),x(%d),%e)',n1,n2,n3,Mosfet(i).Value);
                p{n2,1} = strcat(p{n2,1} , str);
                str = sprintf(' + mosfet_p(x(%d),x(%d),x(%d),%e)',n1,n2,n3,Mosfet(i).Value);
                p{n3,1} = strcat(p{n3,1} , str);
            end

            if (n2==0)
                str = sprintf(' + mosfet_p(x(%d),0,x(%d),%e)',n1,n3,Mosfet(i).Value);
                p{n3,1} = strcat(p{n3,1} , str);
            end
            if (n3==0)
                str = sprintf(' - mosfet_p(x(%d),x(%d),0,%e)',n1,n2,Mosfet(i).Value);
                p{n2,1} = strcat(p{n2,1} , str);
            end 
        end

    end

    for i = 1:sizeX
        p{i} = str2func(p{i});
    end
    func_p = p;
%     B = strcat(regexprep(cellfun(@func2str, p, 'uni', 0), '^@\(x\)', ''), ';');
%     func_p = str2func(strcat('@(x) [', B{:}, ']'));
    
end

function [Nodes, numNode] = addNode(Nodes,Node,numNode)
    try
        Nodes(Node);
    catch
        Nodes(Node) = numNode;
        numNode = numNode + 1;
    end
end


% 
%     % B = zeros(sizeX);
%     % B(1:numNode,1:numNode) = A_C;
%     % 
%     % F = @(x_p,x,t) A*x + E*x_p + eval_fun_array(p,x) + eval_fun_array(r,t);
% %     F = @(x_p,x,t) A*x + E*x_p + func_p(x) + eval_fun_array(r,t);
%     F = @(x_p,x,t) A*x + E*x_p + func_p(x) + func_r(t);
%     % 
%     f2 = @(x) -A*x;
% 
%     f = @(t,x) -A*x - eval_fun_array(p,x) - eval_fun_array(r,t);
% 
%     % x0 = zeros(sizeX,1);
% 
% 
% 
% 
%     % x0 = [4.9996; 4.5941; -0.0046; 0.0000];
% 
%     % x0(4) = 5;
%     % x0(6) = 5;
%     % x0(2) = 0.06183;
%     x0(3) = 5;
%     x0(4) = 0.06183;
%     x0(5) = 5;
%     x0(6) = 5;
%     x0(7) = 0.06183;
%     x0(8) = 0.06183;
%     x0(9) = 5;
%     x0(end-1) = -1.481;
% 
% 
%     t_end =  0.04;
%     N = 1000;
%     h = t_end/N;
%     tol = 1e-8;
%     x_l = x0;
%     x_sol = x0; t_0 = 0; t_sol = t_0;
% 
% 
% 
% 
% 
%     %% Perform dt benchmark
%     % h = h/10;
%     % Joptions.diffvar = 1;
%     % Joptions.vectvars = [];
%     % Joptions.thresh = 1e-3*ones(size(x0));
%     % Joptions.fac = [];
%     % f0 = eval_fun_array(p,x0);        
%     % 
%     % [J_p,Joptions.fac,nF] = odenumjac(@(x) eval_fun_array(p,x), {x0}, f0, Joptions); 
%     % J = A+1/h*E+J_p;
%     % 
%     % [x_bench, error] = NR(F,x0,x0,h,h,J,tol);
%     % 
%     % dt = (x_bench-x0)/h;
%     % dt = dt./(max(dt));
%     % 
%     % h = h/2;
%     % h = t_end/N;
% 
% 
%     % x0(3) = 5;
%     % x0(5) = 5;
%     % x0(4) = 0;
%     % x0(7) = 5;
% 
%     tic;
%     for i = 1: N
%         t_l = t_0 + i*h;
%         l = 1; flag = 1;
% 
%         if i == 1
%             Joptions.diffvar = 1;
%             Joptions.vectvars = [];
%             Joptions.thresh = 1e-3*ones(size(x_l));
%             Joptions.fac = [];
% 
% 
%             f0 = eval_fun_array(p,x_l);        
% 
%             [J_p,Joptions.fac,nF] = odenumjac(@(x) eval_fun_array(p,x), {x_l}, f0, Joptions); 
%             J = A+1/h*E+J_p;
%         end
% 
%         [x_l_1, error] = NR(F,x_l,x0,h,t_l,J,tol);
% 
%         if error > 1e-8 || isnan(error)
%             f0 = eval_fun_array(p,x_l);
%             [J_p,Joptions.fac,nF] = odenumjac(@(x) eval_fun_array(p,x), {x_l}, f0, Joptions); 
%             J = A+1/h*E+J_p;
%             [x_l_1, error] = NR(F,x_l,x0,h,t_l,J,tol);
%         end
%         x_l = x_l_1; x0 = x_l;
%         x_sol = [x_sol x_l];
%         t_sol = [t_sol t_l];
%     end
%     toc
% 
% 
% 
% 
%     figure(); hold on;
%     multiplier = ones(sizeX,1);% multiplier(5) = 1; multiplier(6) = 1;
%     % plot(1000*t_sol,multiplier.*x_sol(4:end,:),'LineWidth',2);
%     plot(1000*t_sol,x_sol(1:end,:),'LineWidth',2);
%     grid on;
%     xlabel('Time in ms');
%     ylabel('Nodal voltages in V and currents in mA');
%     legend
%     % legend('u_1','u_2','u_3','u_4','I_L','I_{V}');
%     legend(keys(Nodes));
%     set(gca, 'FontName', 'Times New Roman','FontSize',14)
% end
% 
% 
% function out = eval_fun_array(fun_ar,x)
%     n=numel(fun_ar);
%     out = zeros(n,1);
%     for i = 1:n
%         out(i) = fun_ar{i}(x);
%     end
% end
% 
% 
% function out = diode(I_s,v1,v2)
%       out = I_s*(exp((v1-v2)/0.0256)-1);
% end
% 
% 
% function g = f_gappy(t,y,mor_object,A,p,r)
%     n = size(y,1);
%     out = zeros(mor_object.g,1);
%     for j = 1:mor_object.g
%         i = mor_object.S(j);
%         out(j,:) = p{i}(y);
%     end
%     
%     f = dot(repmat(out,1,mor_object.g),mor_object.U_gap_S).';
%     b = mor_object.M*f;
%     g_tilde = sum(b'.*mor_object.U_gap,2);    
% %     g = zeros(mor_object.dim_L,1);
%     g = zeros(size(y,1),1);
%     g(mor_object.S) = out;
%     g(mor_object.idx) = g_tilde(mor_object.idx);
%     
%     g = -A*y -g - eval_fun_array(r,t);
%     
% end
% 
% function g = p_gappy(t,y,mor_object,p)
%     n = size(y,1);
%     out = zeros(mor_object.g,1);
%     for j = 1:mor_object.g
%         i = mor_object.S(j);
%         out(j,:) = p{i}(y);
%     end
%     
%     f = dot(repmat(out,1,mor_object.g),mor_object.U_gap_S).';
%     b = mor_object.M*f;
%     g_tilde = sum(b'.*mor_object.U_gap,2);    
% %     g = zeros(mor_object.dim_L,1);
%     g = zeros(size(y,1),1);
%     g(mor_object.S) = out;
%     g(mor_object.idx) = g_tilde(mor_object.idx);    
% end
% 
% function [out, error] = NR(F,x_l,x_0,h,t_l,J,tol)
%     flag = 1; l = 1;
%     while flag
%         x_p = (x_l-x_0)./h;
%         F_val = F(x_p,x_l,t_l);
%         dx = J\(-F_val);
%         x_l_1 = x_l + dx;
%         error = norm(x_l_1 - x_l) + norm(F_val);
%         if isnan(error)
%             break;
%         end
%         
%         if error > tol && l < 30
%             x_l = x_l_1;
%             l = l+1;   
%         else
%             flag = 0;
%         end
%     end
%     
%     out = x_l;
% end
% 
% 

