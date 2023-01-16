function N = nonlinearity(obj)
% assembles nonlinearity in a tensor format
% not optimized

vals00 = zeros(obj.N^3,1);
vals01 = zeros(obj.N^3,1);
subs = zeros(obj.N^3,3);

for p = 1:obj.N
    for j = 1:obj.N
        for m = 1:obj.N
            index = (p-1)*obj.N^2 + (j-1)*obj.N + m;
            subs(index,:) = [p j m];
            vals00(index) = obj.D0(p,j)*obj.D0(p,m);
            vals01(index) = obj.D0(p,j)*obj.D1(p,m);
        end
    end
end

numl = obj.N^3;
bar = obj.bar;
dim = bar*obj.Q;

if obj.TW
    mcell = cell(1,obj.Q);
    for j = 0:obj.Q-1
        mcell{j+1} = -1i*j*obj.B;
    end
    M = toreal(blkdiag(mcell{:}));
    M = sptensor(blkdiag(M,zeros(2)));
    SUBS = zeros(sum(1:(obj.Q-1))*36*4*numl+(obj.Q^2-sum(1:(obj.Q-1)))*18*4*numl+length(M.vals),3); 
    VALS = zeros(sum(1:(obj.Q-1))*36*4*numl+(obj.Q^2-sum(1:(obj.Q-1)))*18*4*numl+length(M.vals),1);
else
    SUBS = zeros(sum(1:(obj.Q-1))*36*4*numl+(obj.Q^2-sum(1:(obj.Q-1)))*18*4*numl,3); 
    VALS = zeros(sum(1:(obj.Q-1))*36*4*numl+(obj.Q^2-sum(1:(obj.Q-1)))*18*4*numl,1);
end
entries_done = 0;
for q = 0:obj.Q-1
for l = 0:obj.Q-1
    n = q-l;
    if n < 0
        interval = entries_done+(1:36*4*numl);
        interval = reshape(interval,[length(interval)/36 36])';

        [SUBS(interval(1,:),:),VALS(interval(1,:))]=...
            sv(dim,subs,obj.Re*n*obj.k*vals00,[q*bar l*bar -n*bar],'im','conj2'); 
        [SUBS(interval(2,:),:),VALS(interval(2,:))]=...
            sv(dim,subs,obj.Re*n*obj.k*vals00,[q*bar+obj.N l*bar -n*bar+obj.N],'im','conj2');
        [SUBS(interval(3,:),:),VALS(interval(3,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar l*bar+obj.N -n*bar],'re','conj2');
        [SUBS(interval(4,:),:),VALS(interval(4,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar+obj.N l*bar+obj.N -n*bar+obj.N],'re','conj2');
        
        [SUBS(interval(5,:),:),VALS(interval(5,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+3*obj.N l*bar -n*bar+3*obj.N],'im','conj2');
        [SUBS(interval(6,:),:),VALS(interval(6,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+4*obj.N l*bar -n*bar+4*obj.N],'im','conj2');
        [SUBS(interval(7,:),:),VALS(interval(7,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+5*obj.N l*bar -n*bar+5*obj.N],'im','conj2');
        
        [SUBS(interval(8,:),:),VALS(interval(8,:))]=...
            sv(dim,subs,vals01,[q*bar+3*obj.N l*bar+obj.N -n*bar+3*obj.N],'re','conj2');
        [SUBS(interval(9,:),:),VALS(interval(9,:))]=...
            sv(dim,subs,vals01,[q*bar+4*obj.N l*bar+obj.N -n*bar+4*obj.N],'re','conj2');
        [SUBS(interval(10,:),:),VALS(interval(10,:))]=...
            sv(dim,subs,vals01,[q*bar+5*obj.N l*bar+obj.N -n*bar+5*obj.N],'re','conj2');
        
        [SUBS(interval(11,:),:),VALS(interval(11,:))]=... 
            sv(dim,subs,-2*n*obj.k*vals00,[q*bar+3*obj.N l*bar+3*obj.N -n*bar],'im','conj2');
        [SUBS(interval(12,:),:),VALS(interval(12,:))]=...
            sv(dim,subs,-2*n*obj.k*vals00,[q*bar+4*obj.N l*bar+5*obj.N -n*bar+obj.N],'im','conj2');
        [SUBS(interval(13,:),:),VALS(interval(13,:))]=...
            sv(dim,subs,-n*obj.k*vals00,[q*bar+5*obj.N l*bar+3*obj.N -n*bar+obj.N],'im','conj2');
        [SUBS(interval(14,:),:),VALS(interval(14,:))]=...
            sv(dim,subs,-n*obj.k*vals00,[q*bar+5*obj.N l*bar+5*obj.N -n*bar],'im','conj2');
      
        [SUBS(interval(15,:),:),VALS(interval(15,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+3*obj.N l*bar+5*obj.N -n*bar],'re','conj2');
        [SUBS(interval(16,:),:),VALS(interval(16,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+4*obj.N l*bar+4*obj.N -n*bar+obj.N],'re','conj2');
        [SUBS(interval(17,:),:),VALS(interval(17,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N l*bar+5*obj.N -n*bar+obj.N],'re','conj2');
        [SUBS(interval(18,:),:),VALS(interval(18,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N l*bar+4*obj.N -n*bar],'re','conj2');
        
        [SUBS(interval(19,:),:),VALS(interval(19,:))]=...
            sv(dim,subs,obj.Re*l*obj.k*vals00,[q*bar -n*bar l*bar],'im','conj1');
        [SUBS(interval(20,:),:),VALS(interval(20,:))]=...
            sv(dim,subs,obj.Re*l*obj.k*vals00,[q*bar+obj.N -n*bar l*bar+obj.N],'im','conj1');
        [SUBS(interval(21,:),:),VALS(interval(21,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar -n*bar+obj.N l*bar],'re','conj1');
        [SUBS(interval(22,:),:),VALS(interval(22,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar+obj.N -n*bar+obj.N l*bar+obj.N],'re','conj1');
        
        [SUBS(interval(23,:),:),VALS(interval(23,:))]=...
            sv(dim,subs,l*obj.k*vals00,[q*bar+3*obj.N -n*bar l*bar+3*obj.N],'im','conj1');
        [SUBS(interval(24,:),:),VALS(interval(24,:))]=...
            sv(dim,subs,l*obj.k*vals00,[q*bar+4*obj.N -n*bar l*bar+4*obj.N],'im','conj1');
        [SUBS(interval(25,:),:),VALS(interval(25,:))]=...
            sv(dim,subs,l*obj.k*vals00,[q*bar+5*obj.N -n*bar l*bar+5*obj.N],'im','conj1');
        
        [SUBS(interval(26,:),:),VALS(interval(26,:))]=...
            sv(dim,subs,vals01,[q*bar+3*obj.N -n*bar+obj.N l*bar+3*obj.N],'re','conj1');
        [SUBS(interval(27,:),:),VALS(interval(27,:))]=...
            sv(dim,subs,vals01,[q*bar+4*obj.N -n*bar+obj.N l*bar+4*obj.N],'re','conj1');
        [SUBS(interval(28,:),:),VALS(interval(28,:))]=...
            sv(dim,subs,vals01,[q*bar+5*obj.N -n*bar+obj.N l*bar+5*obj.N],'re','conj1');
        
        [SUBS(interval(29,:),:),VALS(interval(29,:))]=... 
            sv(dim,subs,-2*l*obj.k*vals00,[q*bar+3*obj.N -n*bar+3*obj.N l*bar],'im','conj1');
        [SUBS(interval(30,:),:),VALS(interval(30,:))]=...
            sv(dim,subs,-2*l*obj.k*vals00,[q*bar+4*obj.N -n*bar+5*obj.N l*bar+obj.N],'im','conj1');
        [SUBS(interval(31,:),:),VALS(interval(31,:))]=...
            sv(dim,subs,-l*obj.k*vals00,[q*bar+5*obj.N -n*bar+3*obj.N l*bar+obj.N],'im','conj1');
        [SUBS(interval(32,:),:),VALS(interval(32,:))]=...
            sv(dim,subs,-l*obj.k*vals00,[q*bar+5*obj.N -n*bar+5*obj.N l*bar],'im','conj1');
      
        [SUBS(interval(33,:),:),VALS(interval(33,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+3*obj.N -n*bar+5*obj.N l*bar],'re','conj1');
        [SUBS(interval(34,:),:),VALS(interval(34,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+4*obj.N -n*bar+4*obj.N l*bar+obj.N],'re','conj1');
        [SUBS(interval(35,:),:),VALS(interval(35,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N -n*bar+5*obj.N l*bar+obj.N],'re','conj1');
        [SUBS(interval(36,:),:),VALS(interval(36,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N -n*bar+4*obj.N l*bar],'re','conj1');
        
        entries_done = entries_done + 36*4*numl;
    else
        interval = entries_done+(1:18*4*numl);
        interval = reshape(interval,[length(interval)/18 18])';
        
        [SUBS(interval(1,:),:),VALS(interval(1,:))]=...
            sv(dim,subs,obj.Re*n*obj.k*vals00,[q*bar l*bar n*bar],'im','std');
        [SUBS(interval(2,:),:),VALS(interval(2,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar l*bar+obj.N n*bar],'re','std');
        [SUBS(interval(3,:),:),VALS(interval(3,:))]=...
            sv(dim,subs,obj.Re*n*obj.k*vals00,[q*bar+obj.N l*bar n*bar+obj.N],'im','std');
        [SUBS(interval(4,:),:),VALS(interval(4,:))]=...
            sv(dim,subs,obj.Re*vals01,[q*bar+obj.N l*bar+obj.N n*bar+obj.N],'re','std');
        
        [SUBS(interval(5,:),:),VALS(interval(5,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+3*obj.N l*bar n*bar+3*obj.N],'im','std');
        [SUBS(interval(6,:),:),VALS(interval(6,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+4*obj.N l*bar n*bar+4*obj.N],'im','std');
        [SUBS(interval(7,:),:),VALS(interval(7,:))]=...
            sv(dim,subs,n*obj.k*vals00,[q*bar+5*obj.N l*bar n*bar+5*obj.N],'im','std');
        
        [SUBS(interval(8,:),:),VALS(interval(8,:))]=...
            sv(dim,subs,vals01,[q*bar+3*obj.N l*bar+obj.N n*bar+3*obj.N],'re','std');
        [SUBS(interval(9,:),:),VALS(interval(9,:))]=...
            sv(dim,subs,vals01,[q*bar+4*obj.N l*bar+obj.N n*bar+4*obj.N],'re','std');
        [SUBS(interval(10,:),:),VALS(interval(10,:))]=...
            sv(dim,subs,vals01,[q*bar+5*obj.N l*bar+obj.N n*bar+5*obj.N],'re','std');
        
        [SUBS(interval(11,:),:),VALS(interval(11,:))]=... 
            sv(dim,subs,-2*n*obj.k*vals00,[q*bar+3*obj.N l*bar+3*obj.N n*bar],'im','std');
        [SUBS(interval(12,:),:),VALS(interval(12,:))]=...
            sv(dim,subs,-2*n*obj.k*vals00,[q*bar+4*obj.N l*bar+5*obj.N n*bar+obj.N],'im','std');
        [SUBS(interval(13,:),:),VALS(interval(13,:))]=...
            sv(dim,subs,-n*obj.k*vals00,[q*bar+5*obj.N l*bar+3*obj.N n*bar+obj.N],'im','std');
        [SUBS(interval(14,:),:),VALS(interval(14,:))]=...
            sv(dim,subs,-n*obj.k*vals00,[q*bar+5*obj.N l*bar+5*obj.N n*bar],'im','std');
      
        [SUBS(interval(15,:),:),VALS(interval(15,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+3*obj.N l*bar+5*obj.N n*bar],'re','std');
        [SUBS(interval(16,:),:),VALS(interval(16,:))]=...
            sv(dim,subs,-2*vals01,[q*bar+4*obj.N l*bar+4*obj.N n*bar+obj.N],'re','std');
        [SUBS(interval(17,:),:),VALS(interval(17,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N l*bar+5*obj.N n*bar+obj.N],'re','std');
        [SUBS(interval(18,:),:),VALS(interval(18,:))]=...
            sv(dim,subs,-vals01,[q*bar+5*obj.N l*bar+4*obj.N n*bar],'re','std');
 
        entries_done = entries_done + 18*4*numl;
    end
end
end

if obj.TW
    SUBS(entries_done+1:end,:) = [M.subs obj.dim*ones(length(M.vals),1)];
    VALS(entries_done+1:end) = M.vals;
end

N = sptensor(SUBS,VALS,[obj.dim obj.dim obj.dim]);
N = obj.stopinterf(N);
N = obj.remove_extra_lines(N);
end

function [sabs,vels] = sv(dim,subs,vals,stpos,imre,spec)
% RexRe
% RexIm
% ImxRe
% ImxIm
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
end
numl = length(vals);
vels = zeros(4*numl,1);
sabs = zeros(4*numl,3);
switch imre
    case 'im'
        sabs(1:numl,:) = repmat([stpos(1)+dim stpos(2:3)],[numl 1])+subs; 
        vels(1:numl) = con(1)*vals; 
        sabs(numl+1:2*numl,:) = repmat([stpos(1:2) stpos(3)+dim],[numl 1])+subs; 
        vels(numl+1:2*numl) = -con(2)*vals; 
        sabs(2*numl+1:3*numl,:) = repmat([stpos(1) stpos(2)+dim stpos(3)],[numl 1])+subs; 
        vels(2*numl+1:3*numl) = -con(3)*vals; 
        sabs(3*numl+1:end,:) = repmat(stpos+dim,[numl 1])+subs; 
        vels(3*numl+1:end) = -con(4)*vals; 
    case 're'
        sabs(1:numl,:) = repmat(stpos,[numl 1])+subs; 
        vels(1:numl) = con(1)*vals; 
        sabs(numl+1:2*numl,:) = repmat([stpos(1)+dim stpos(2) stpos(3)+dim],[numl 1])+subs; 
        vels(numl+1:2*numl) = con(2)*vals; 
        sabs(2*numl+1:3*numl,:) = repmat([stpos(1:2)+dim stpos(3)],[numl 1])+subs; 
        vels(2*numl+1:3*numl) = con(3)*vals; 
        sabs(3*numl+1:end,:) = repmat([stpos(1) stpos(2:3)+dim],[numl 1])+subs; 
        vels(3*numl+1:end) = -con(4)*vals; 
end
end