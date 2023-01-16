function N = nonlinearity_multiindex(obj)
% assembles nonlinearity in a multiindex format (fast)
% Below lie the worst lines of spaghetti you've ever witnessed.
% For the sake of your own sanity, do not read.

% numl = obj.N^2*(obj.N-2);
numl = obj.N^2*(obj.N);
vals00 = zeros(numl,1);
vals01 = zeros(numl,1);
vals00t = zeros(numl,1);
vals01t = zeros(numl,1);
vals_subs = zeros(numl,2);
mults_init = zeros(obj.N^2,2); 
numle = obj.N^2;

vals_subsh = zeros(sum(1:obj.N)*(obj.N-2),2);
valsh00 = zeros(sum(1:obj.N)*(obj.N-2),1);
valsh01 = zeros(sum(1:obj.N)*(obj.N-2),1);
multsh = zeros(sum(1:obj.N),2);
for p = 1:obj.N 
    for j = 1:obj.N
        for m = 1:obj.N
            index = (p-1)*obj.N^2 + (j-1)*obj.N + m;
            index2 = (j-1)*obj.N + m;
            mults_init(index2,:) = [j m];
            vals_subs(index,:) = [p index2];
            vals00(index) = obj.D0(p,j)*obj.D0(p,m);
            vals00t(index) = obj.D0(p,m)*obj.D0(p,j);
            vals01(index) = obj.D0(p,j)*obj.D1(p,m);
            vals01t(index) = obj.D0(p,m)*obj.D1(p,j);
        end
        for m = j:obj.N
            index = (p-1)*sum(1:obj.N) + sum((obj.N-j+2):obj.N) + m+1-j;
            index2 = sum((obj.N-j+2):obj.N) + m+1-j; 
            multsh(index2,:) = [j m];
            vals_subsh(index,:) = [p index2];
            if j == m
            valsh00(index) = obj.D0(p,j)*obj.D0(p,m);
            valsh01(index) = obj.D0(p,j)*obj.D1(p,m);
            else
            valsh00(index) = obj.D0(p,j)*obj.D0(p,m) + obj.D0(p,m)*obj.D0(p,j);
            valsh01(index) = obj.D0(p,j)*obj.D1(p,m) + obj.D0(p,m)*obj.D1(p,j);
            end
        end
    end
end

bar = obj.bar;
dim = bar*obj.Q-bar;
k = obj.k;

if obj.TW
mcell = cell(1,obj.Q);
B = obj.stopinterf(obj.B);
for j = 0:obj.Q-1
    mcell{j+1} = -1i*j*B;
end
M = toreal(blkdiag(mcell{:}));
M = sptensor(blkdiag(M,zeros(2)));
M = obj.remove_extra_lines(M);
[sabs,is] = sort(M.subs(:,2));
isc = diff([0; sabs])~=0;
MULTS = zeros(sum(isc)+16*4*numle*(sum(1:obj.Q-2))+4*numle*12*(obj.Q-1)+2*numle*11*(obj.Q-1)+4*numle*ceil((obj.Q-1)/2)+4*numle*2*floor((obj.Q-1)/2)+4*numle*8*floor((obj.Q-1)/2)+2*1*numle,2);
VALSUBS = zeros(length(M.vals)+26*4*numl*(sum(1:obj.Q-2))+26*4*numl*(sum((obj.Q-1):-2:1)-(obj.Q-1))+4*numl*21*(obj.Q-1)+2*numl*17*(obj.Q-1)+4*numl*2*floor((obj.Q-1)/2)+4*numl*12*floor((obj.Q-1)/2)+2*numl,2);
VAL = zeros(length(M.vals)+26*4*numl*(sum(1:obj.Q-2))+26*4*numl*(sum((obj.Q-1):-2:1)-(obj.Q-1))+4*numl*21*(obj.Q-1)+2*numl*17*(obj.Q-1)+4*numl*2*floor((obj.Q-1)/2)+4*numl*12*floor((obj.Q-1)/2)+2*numl,1);
else
MULTS = zeros(16*4*numle*(sum(1:obj.Q-2))+4*numle*12*(obj.Q-1)+2*numle*11*(obj.Q-1)+4*numle*ceil((obj.Q-1)/2)+4*numle*2*floor((obj.Q-1)/2)+4*numle*8*floor((obj.Q-1)/2)+2*1*numle,2);
VALSUBS = zeros(26*4*numl*(sum(1:obj.Q-2))+26*4*numl*(sum((obj.Q-1):-2:1)-(obj.Q-1))+4*numl*21*(obj.Q-1)+2*numl*17*(obj.Q-1)+4*numl*2*floor((obj.Q-1)/2)+4*numl*12*floor((obj.Q-1)/2)+2*numl,2);
VAL = zeros(26*4*numl*(sum(1:obj.Q-2))+26*4*numl*(sum((obj.Q-1):-2:1)-(obj.Q-1))+4*numl*21*(obj.Q-1)+2*numl*17*(obj.Q-1)+4*numl*2*floor((obj.Q-1)/2)+4*numl*12*floor((obj.Q-1)/2)+2*numl,1);
end

MIentries_done = 0;
Ventries_done = 0;

n = 0;
l = 0;
q = 0;
Mints = MIentries_done+(1:2*1*numle); 
Mints = reshape(Mints,[length(Mints)/2 2])';
Vints = Ventries_done+(1:2*1*numl); 
Vints = reshape(Vints,[length(Vints)/2 2])';
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser_nq(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar+2*obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 1*numle;
[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser_nq(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar+2*obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 1*numle;
Ventries_done = Ventries_done + 2*1*numl;


for l = 1:obj.Q-1
q = l; 
Mints = MIentries_done+(1:11*2*numle); 
Mints = reshape(Mints,[length(Mints)/11 11])';
Vints = Ventries_done+(1:17*2*numl); 
Vints = reshape(Vints,[length(Vints)/17 17])';

[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser_n(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar n*bar+2*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;

[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser_n(dim,obj.Re*vals01t,obj.Re*vals01,[q*bar n*bar+2*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(3,:)),VALSUBS(Vints(3,:),:),~] =... 
vels_parser_n(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar+obj.N n*bar+2*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;

Ventries_done = Ventries_done + 3*2*numl;

%q3 1 + 3 
[VAL(Vints(4,:)),VALSUBS(Vints(4,:),:),MULTS(Mints(3,:),:)] =... 
vels_parser_n(dim,l*k*vals00,l*k*vals00t,[q*bar+3*obj.N n*bar+2*obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
[VAL(Vints(5,:)),VALSUBS(Vints(5,:),:),MULTS(Mints(4,:),:)] =... 
vels_parser_n(dim,-2*l*k*vals00,-2*l*k*vals00t,[q*bar+3*obj.N n*bar+3*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 2*2*numl;
% q3 2 + q6 3
[VAL(Vints(6,:)),VALSUBS(Vints(6,:),:),MULTS(Mints(5,:),:)] =... 
vels_parser_n(dim,vals01t,vals01,[q*bar+3*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(7,:)),VALSUBS(Vints(7,:),:),~] =... 
vels_parser_n(dim,-l*k*vals00,-l*k*vals00t,[q*bar+5*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 2*2*numl;

% 4, also 1 and 4 of q6
[VAL(Vints(8,:)),VALSUBS(Vints(8,:),:),MULTS(Mints(6,:),:)] =... 
vels_parser_n(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar+2*obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(9,:)),VALSUBS(Vints(9,:),:),~] =... 
vels_parser_n(dim,l*k*vals00,l*k*vals00t,[q*bar+5*obj.N n*bar+2*obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
[VAL(Vints(10,:)),VALSUBS(Vints(10,:),:),MULTS(Mints(7,:),:)] =... 
vels_parser_n(dim,-2*vals01,-2*vals01t,[q*bar+3*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(11,:)),VALSUBS(Vints(11,:),:),~] =... 
vels_parser_n(dim,-l*k*vals00t,-l*k*vals00,[q*bar+5*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 4*2*numl;

%q5
% 1 + q6 6
[VAL(Vints(12,:)),VALSUBS(Vints(12,:),:),MULTS(Mints(8,:),:)] =... 
vels_parser_n(dim,l*k*vals00,l*k*vals00t,[q*bar+4*obj.N n*bar+2*obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(13,:)),VALSUBS(Vints(13,:),:),~] =... 
vels_parser_n(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar+2*obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 2*numle;
[VAL(Vints(14,:)),VALSUBS(Vints(14,:),:),MULTS(Mints(9,:),:)] =... 
vels_parser_n(dim,-vals01,-vals01t,[q*bar+5*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 3*2*numl;

%q5 2 + 4
[VAL(Vints(15,:)),VALSUBS(Vints(15,:),:),MULTS(Mints(10,:),:)] =... 
vels_parser_n(dim,vals01t-2*vals01,vals01-2*vals01t,[q*bar+4*obj.N n*bar+4*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 1*2*numl;

%q5 3 + q6 2 5
[VAL(Vints(16,:)),VALSUBS(Vints(16,:),:),MULTS(Mints(11,:),:)] =... 
vels_parser_n(dim,vals01t-vals01,vals01-vals01t,[q*bar+5*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(17,:)),VALSUBS(Vints(17,:),:),~] =... 
vels_parser_n(dim,-2*l*k*vals00,-2*l*k*vals00t,[q*bar+4*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 2*numle;
Ventries_done = Ventries_done + 2*2*numl;


end

for n = 1:obj.Q-1
l = n;
q = l+n;
if q < obj.Q
Mints = MIentries_done+(1:2*(4*numle));
Mints = reshape(Mints,[length(Mints)/2 2])';
Vints = Ventries_done+(1:2*(4*numl));
Vints = reshape(Vints,[length(Vints)/2 2])';
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar n*bar l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser(dim,obj.Re*vals01,obj.Re*vals01t,[q*bar n*bar l*bar]+obj.N,...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

Mints = MIentries_done+(1:8*4*numle);
Mints = reshape(Mints,[length(Mints)/8 8])';
Vints = Ventries_done+(1:12*4*numl);
Vints = reshape(Vints,[length(Vints)/12 12])';

[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar+obj.N n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser(dim,obj.Re*vals01,obj.Re*vals01t,[q*bar  n*bar+obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;


%q3 1 + 3 
[VAL(Vints(3,:)),VALSUBS(Vints(3,:),:),MULTS(Mints(3,:),:)] =... 
vels_parser(dim,l*k*vals00-2*n*k*vals00t,l*k*vals00t-2*n*k*vals00,[q*bar+3*obj.N n*bar l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 1*4*numl;
% q3 2 + q6 3
[VAL(Vints(4,:)),VALSUBS(Vints(4,:),:),MULTS(Mints(4,:),:)] =... 
vels_parser(dim,vals01t,vals01,[q*bar+3*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(5,:)),VALSUBS(Vints(5,:),:),~] =... 
vels_parser(dim,-l*k*vals00,-l*k*vals00t,[q*bar+5*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

% 4, also 1 and 4 of q6
[VAL(Vints(6,:)),VALSUBS(Vints(6,:),:),MULTS(Mints(5,:),:)] =... 
vels_parser(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(7,:)),VALSUBS(Vints(7,:),:),~] =... 
vels_parser(dim,l*k*vals00-n*k*vals00t,l*k*vals00t-n*k*vals00,[q*bar+5*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

%q5
% 1 + q6 6
[VAL(Vints(8,:)),VALSUBS(Vints(8,:),:),MULTS(Mints(6,:),:)] =... 
vels_parser(dim,l*k*vals00,l*k*vals00t,[q*bar+4*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(9,:)),VALSUBS(Vints(9,:),:),~] =... 
vels_parser(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

%q5 2 + 4
[VAL(Vints(10,:)),VALSUBS(Vints(10,:),:),MULTS(Mints(7,:),:)] =... 
vels_parser(dim,vals01-2*vals01t,vals01t-2*vals01,[q*bar+4*obj.N n*bar+obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 1*4*numl;

%q5 3 + q6 2 5
[VAL(Vints(11,:)),VALSUBS(Vints(11,:),:),MULTS(Mints(8,:),:)] =... 
vels_parser(dim,-2*n*k*vals00t,-2*n*k*vals00,[q*bar+4*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(12,:)),VALSUBS(Vints(12,:),:),~] =... 
vels_parser(dim,vals01-vals01t,vals01t-vals01,[q*bar+5*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;


q = 0;
MIentries_done = MIentries_done - 8*4*numle;
Vints = Ventries_done+(1:1*4*numl);
Vints = reshape(Vints,[length(Vints)/1 1])';
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),~] =... 
vels_parser_q(dim,obj.Re*vals01t*2,obj.Re*vals01*2,[q*bar+2*obj.N n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 8*4*numle;
Ventries_done = Ventries_done + 1*4*numl;

else
Mints = MIentries_done+(1:1*4*numle);
Mints = reshape(Mints,[length(Mints)/1 1])';
Vints = Ventries_done+(1:1*4*numl); 
Vints = reshape(Vints,[length(Vints)/1 1])';
q = 0;
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser_q(dim,obj.Re*vals01t*2,obj.Re*vals01*2,[q*bar+2*obj.N n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1'); 
% the indices here are actually conj so it zeros out where it needs to
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 1*4*numl;
end

Mints = MIentries_done+(1:12*4*numle);
Mints = reshape(Mints,[length(Mints)/12 12])';
Vints = Ventries_done+(1:20*4*numl); 
Vints = reshape(Vints,[length(Vints)/20 20])';

%q3 1 + 3 
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =... 
vels_parser_q(dim,l*k*vals00+2*n*k*vals00t,l*k*vals00t+2*n*k*vals00,[q*bar+3*obj.N n*bar l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser_q(dim,-n*k*vals00t-2*l*k*vals00,-n*k*vals00-2*l*k*vals00t,[q*bar+3*obj.N n*bar+3*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;
% q3 2 + q6 3
[VAL(Vints(3,:)),VALSUBS(Vints(3,:),:),MULTS(Mints(3,:),:)] =... 
vels_parser_q(dim,vals01t,vals01,[q*bar+3*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(4,:)),VALSUBS(Vints(4,:),:),~] =... 
vels_parser_q(dim,-l*k*vals00,-l*k*vals00t,[q*bar+5*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(5,:)),VALSUBS(Vints(5,:),:),MULTS(Mints(4,:),:)] =... 
vels_parser_q(dim,vals01,vals01t,[q*bar+3*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(6,:)),VALSUBS(Vints(6,:),:),~] =... 
vels_parser_q(dim,n*k*vals00t,n*k*vals00,[q*bar+5*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

% 4, also 1 and 4 of q6
[VAL(Vints(7,:)),VALSUBS(Vints(7,:),:),MULTS(Mints(5,:),:)] =... 
vels_parser_q(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(8,:)),VALSUBS(Vints(8,:),:),~] =... 
vels_parser_q(dim,l*k*vals00+n*k*vals00t,l*k*vals00t+n*k*vals00,[q*bar+5*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(9,:)),VALSUBS(Vints(9,:),:),MULTS(Mints(6,:),:)] =... 
vels_parser_q(dim,-2*vals01,-2*vals01t,[q*bar+3*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(10,:)),VALSUBS(Vints(10,:),:),~] =... 
vels_parser_q(dim,-n*k*vals00-l*k*vals00t,-n*k*vals00t-l*k*vals00,[q*bar+5*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

%q5
% 1 + q6 6
[VAL(Vints(11,:)),VALSUBS(Vints(11,:),:),MULTS(Mints(7,:),:)] =... 
vels_parser_q(dim,l*k*vals00,l*k*vals00t,[q*bar+4*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(12,:)),VALSUBS(Vints(12,:),:),~] =... 
vels_parser_q(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(13,:)),VALSUBS(Vints(13,:),:),MULTS(Mints(8,:),:)] =... 
vels_parser_q(dim,-n*k*vals00t,-n*k*vals00,[q*bar+4*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(14,:)),VALSUBS(Vints(14,:),:),~] =... 
vels_parser_q(dim,-vals01,-vals01t,[q*bar+5*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

%q5 2 + 4
[VAL(Vints(15,:)),VALSUBS(Vints(15,:),:),MULTS(Mints(9,:),:)] =... 
vels_parser_q(dim,vals01-2*vals01t,vals01t-2*vals01,[q*bar+4*obj.N n*bar+obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(16,:)),VALSUBS(Vints(16,:),:),MULTS(Mints(10,:),:)] =... 
vels_parser_q(dim,vals01t-2*vals01,vals01-2*vals01t,[q*bar+4*obj.N n*bar+4*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

%q5 3 + q6 2 5
[VAL(Vints(17,:)),VALSUBS(Vints(17,:),:),MULTS(Mints(11,:),:)] =... 
vels_parser_q(dim,2*n*k*vals00t,2*n*k*vals00,[q*bar+4*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(18,:)),VALSUBS(Vints(18,:),:),~] =... 
vels_parser_q(dim,vals01-vals01t,vals01t-vals01,[q*bar+5*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(19,:)),VALSUBS(Vints(19,:),:),MULTS(Mints(12,:),:)] =... 
vels_parser_q(dim,-2*l*k*vals00,-2*l*k*vals00t,[q*bar+4*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(20,:)),VALSUBS(Vints(20,:),:),~] =... 
vels_parser_q(dim,vals01t-vals01,vals01-vals01t,[q*bar+5*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;


for l = n+1:obj.Q-1 

Mints = MIentries_done+(1:16*4*numle);
Mints = reshape(Mints,[length(Mints)/16 16])';

if n+l < obj.Q
    Vints = Ventries_done+(1:52*4*numl);
    Vints = reshape(Vints,[length(Vints)/52 52])';
else
    Vints = Ventries_done+(1:26*4*numl);
    Vints = reshape(Vints,[length(Vints)/26 26])';
end

q = l-n;
[VAL(Vints(1,:)),VALSUBS(Vints(1,:),:),MULTS(Mints(1,:),:)] =...
vels_parser(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar n*bar l*bar],...
MIentries_done,vals_subs,mults_init,-obj.Re*n*k*vals00,-obj.Re*n*k*vals00t,'double');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(2,:)),VALSUBS(Vints(2,:),:),MULTS(Mints(2,:),:)] =... 
vels_parser(dim,obj.Re*vals01,obj.Re*vals01t,[q*bar n*bar l*bar]+obj.N,...
MIentries_done,vals_subs,mults_init,'double2','double2');
MIentries_done = MIentries_done + 4*numle;

[VAL(Vints(3,:)),VALSUBS(Vints(3,:),:),MULTS(Mints(3,:),:)] =... 
vels_parser(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar+obj.N n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(4,:)),VALSUBS(Vints(4,:),:),~] =... 
vels_parser(dim,obj.Re*vals01t,obj.Re*vals01,[q*bar n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(5,:)),VALSUBS(Vints(5,:),:),MULTS(Mints(4,:),:)] =... 
vels_parser(dim,-obj.Re*n*k*vals00t,-obj.Re*n*k*vals00,[q*bar+obj.N n*bar+obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(6,:)),VALSUBS(Vints(6,:),:),~] =... 
vels_parser(dim,obj.Re*vals01,obj.Re*vals01t,[q*bar n*bar+obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 6*4*numl;

%q3 1 + 3 
[VAL(Vints(7,:)),VALSUBS(Vints(7,:),:),MULTS(Mints(5,:),:)] =... 
vels_parser(dim,l*k*vals00+2*n*k*vals00t,l*k*vals00t+2*n*k*vals00,[q*bar+3*obj.N n*bar l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(8,:)),VALSUBS(Vints(8,:),:),MULTS(Mints(6,:),:)] =... 
vels_parser(dim,-n*k*vals00t-2*l*k*vals00,-n*k*vals00-2*l*k*vals00t,[q*bar+3*obj.N n*bar+3*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;
% q3 2 + q6 3
[VAL(Vints(9,:)),VALSUBS(Vints(9,:),:),MULTS(Mints(7,:),:)] =... 
vels_parser(dim,vals01t,vals01,[q*bar+3*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(10,:)),VALSUBS(Vints(10,:),:),~] =... 
vels_parser(dim,-l*k*vals00,-l*k*vals00t,[q*bar+5*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(11,:)),VALSUBS(Vints(11,:),:),MULTS(Mints(8,:),:)] =... 
vels_parser(dim,vals01,vals01t,[q*bar+3*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(12,:)),VALSUBS(Vints(12,:),:),~] =... 
vels_parser(dim,n*k*vals00t,n*k*vals00,[q*bar+5*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

% 4, also 1 and 4 of q6
[VAL(Vints(13,:)),VALSUBS(Vints(13,:),:),MULTS(Mints(9,:),:)] =... 
vels_parser(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(14,:)),VALSUBS(Vints(14,:),:),~] =... 
vels_parser(dim,l*k*vals00+n*k*vals00t,l*k*vals00t+n*k*vals00,[q*bar+5*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1'); 
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(15,:)),VALSUBS(Vints(15,:),:),MULTS(Mints(10,:),:)] =... 
vels_parser(dim,-2*vals01,-2*vals01t,[q*bar+3*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','conj1');
[VAL(Vints(16,:)),VALSUBS(Vints(16,:),:),~] =... 
vels_parser(dim,-n*k*vals00-l*k*vals00t,-n*k*vals00t-l*k*vals00,[q*bar+5*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

%q5
% 1 + q6 6
[VAL(Vints(17,:)),VALSUBS(Vints(17,:),:),MULTS(Mints(11,:),:)] =... 
vels_parser(dim,l*k*vals00,l*k*vals00t,[q*bar+4*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(18,:)),VALSUBS(Vints(18,:),:),~] =... 
vels_parser(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(19,:)),VALSUBS(Vints(19,:),:),MULTS(Mints(12,:),:)] =... 
vels_parser(dim,-n*k*vals00t,-n*k*vals00,[q*bar+4*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(20,:)),VALSUBS(Vints(20,:),:),~] =... 
vels_parser(dim,-vals01,-vals01t,[q*bar+5*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

%q5 2 + 4
[VAL(Vints(21,:)),VALSUBS(Vints(21,:),:),MULTS(Mints(13,:),:)] =... 
vels_parser(dim,vals01-2*vals01t,vals01t-2*vals01,[q*bar+4*obj.N n*bar+obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(22,:)),VALSUBS(Vints(22,:),:),MULTS(Mints(14,:),:)] =... 
vels_parser(dim,vals01t-2*vals01,vals01-2*vals01t,[q*bar+4*obj.N n*bar+4*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

%q5 3 + q6 2 5
[VAL(Vints(23,:)),VALSUBS(Vints(23,:),:),MULTS(Mints(15,:),:)] =... 
vels_parser(dim,2*n*k*vals00t,2*n*k*vals00,[q*bar+4*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(24,:)),VALSUBS(Vints(24,:),:),~] =... 
vels_parser(dim,vals01-vals01t,vals01t-vals01,[q*bar+5*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1'); 
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(25,:)),VALSUBS(Vints(25,:),:),MULTS(Mints(16,:),:)] =... 
vels_parser(dim,-2*l*k*vals00,-2*l*k*vals00t,[q*bar+4*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','conj1');
[VAL(Vints(26,:)),VALSUBS(Vints(26,:),:),~] =... 
vels_parser(dim,vals01t-vals01,vals01-vals01t,[q*bar+5*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','conj1'); 
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;


q = n + l;
if q < obj.Q
MIentries_done = MIentries_done - 16*4*numle;
[VAL(Vints(27,:)),VALSUBS(Vints(27,:),:),~] =... 
vels_parser(dim,obj.Re*n*k*vals00t+obj.Re*l*k*vals00,obj.Re*n*k*vals00+obj.Re*l*k*vals00t,[q*bar n*bar l*bar],... 
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(28,:)),VALSUBS(Vints(28,:),:),~] =... 
vels_parser(dim,obj.Re*vals01+obj.Re*vals01t,obj.Re*vals01t+obj.Re*vals01,[q*bar n*bar l*bar]+obj.N,...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;

[VAL(Vints(29,:)),VALSUBS(Vints(29,:),:),~] =... 
vels_parser(dim,obj.Re*vals01t,obj.Re*vals01,[q*bar n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(30,:)),VALSUBS(Vints(30,:),:),~] =... 
vels_parser(dim,obj.Re*l*k*vals00,obj.Re*l*k*vals00t,[q*bar+obj.N n*bar l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;

[VAL(Vints(31,:)),VALSUBS(Vints(31,:),:),~] =... 
vels_parser(dim,obj.Re*vals01,obj.Re*vals01t,[q*bar  n*bar+obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(32,:)),VALSUBS(Vints(32,:),:),~] =... 
vels_parser(dim,obj.Re*n*k*vals00t,obj.Re*n*k*vals00,[q*bar+obj.N n*bar+obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 6*4*numl;

%q3 1 + 3 
[VAL(Vints(33,:)),VALSUBS(Vints(33,:),:),~] =... 
vels_parser(dim,l*k*vals00-2*n*k*vals00t,l*k*vals00t-2*n*k*vals00,[q*bar+3*obj.N n*bar l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(34,:)),VALSUBS(Vints(34,:),:),~] =... 
vels_parser(dim,+n*k*vals00t-2*l*k*vals00,+n*k*vals00-2*l*k*vals00t,[q*bar+3*obj.N n*bar+3*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;
% q3 2 + q6 3
[VAL(Vints(35,:)),VALSUBS(Vints(35,:),:),~] =... 
vels_parser(dim,vals01t,vals01,[q*bar+3*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(36,:)),VALSUBS(Vints(36,:),:),~] =... 
vels_parser(dim,-l*k*vals00,-l*k*vals00t,[q*bar+5*obj.N n*bar+3*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(37,:)),VALSUBS(Vints(37,:),:),~] =... 
vels_parser(dim,vals01,vals01t,[q*bar+3*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(38,:)),VALSUBS(Vints(38,:),:),~] =... 
vels_parser(dim,-n*k*vals00t,-n*k*vals00,[q*bar+5*obj.N n*bar+obj.N l*bar+3*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

% 4, also 1 and 4 of q6
[VAL(Vints(39,:)),VALSUBS(Vints(39,:),:),~] =... 
vels_parser(dim,-2*vals01t,-2*vals01,[q*bar+3*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(40,:)),VALSUBS(Vints(40,:),:),~] =... 
vels_parser(dim,l*k*vals00-n*k*vals00t,l*k*vals00t-n*k*vals00,[q*bar+5*obj.N n*bar l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(41,:)),VALSUBS(Vints(41,:),:),~] =... 
vels_parser(dim,-2*vals01,-2*vals01t,[q*bar+3*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
[VAL(Vints(42,:)),VALSUBS(Vints(42,:),:),~] =... 
vels_parser(dim,+n*k*vals00-l*k*vals00t,+n*k*vals00t-l*k*vals00,[q*bar+5*obj.N n*bar+5*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;


%q5
% 1 + q6 6
[VAL(Vints(43,:)),VALSUBS(Vints(43,:),:),~] =... 
vels_parser(dim,l*k*vals00,l*k*vals00t,[q*bar+4*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(44,:)),VALSUBS(Vints(44,:),:),~] =... 
vels_parser(dim,-vals01t,-vals01,[q*bar+5*obj.N n*bar l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(45,:)),VALSUBS(Vints(45,:),:),~] =... 
vels_parser(dim,n*k*vals00t,n*k*vals00,[q*bar+4*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(46,:)),VALSUBS(Vints(46,:),:),~] =... 
vels_parser(dim,-vals01,-vals01t,[q*bar+5*obj.N n*bar+4*obj.N l*bar],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

%q5 2 + 4
[VAL(Vints(47,:)),VALSUBS(Vints(47,:),:),~] =... 
vels_parser(dim,vals01-2*vals01t,vals01t-2*vals01,[q*bar+4*obj.N n*bar+obj.N l*bar+4*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(48,:)),VALSUBS(Vints(48,:),:),~] =... 
vels_parser(dim,vals01t-2*vals01,vals01-2*vals01t,[q*bar+4*obj.N n*bar+4*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 2*4*numl;

%q5 3 + q6 2 5
[VAL(Vints(49,:)),VALSUBS(Vints(49,:),:),~] =... 
vels_parser(dim,-2*n*k*vals00t,-2*n*k*vals00,[q*bar+4*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');
[VAL(Vints(50,:)),VALSUBS(Vints(50,:),:),~] =... 
vels_parser(dim,vals01-vals01t,vals01t-vals01,[q*bar+5*obj.N n*bar+obj.N l*bar+5*obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
[VAL(Vints(51,:)),VALSUBS(Vints(51,:),:),~] =... 
vels_parser(dim,-2*l*k*vals00,-2*l*k*vals00t,[q*bar+4*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'im','std');    
[VAL(Vints(52,:)),VALSUBS(Vints(52,:),:),~] =... 
vels_parser(dim,vals01t-vals01,vals01-vals01t,[q*bar+5*obj.N n*bar+5*obj.N l*bar+obj.N],...
MIentries_done,vals_subs,mults_init,'re','std');
MIentries_done = MIentries_done + 4*numle;
Ventries_done = Ventries_done + 4*4*numl;

end



end
end

VALSUBS(:,1) = VALSUBS(:,1) - 2*obj.N;
MULTS = MULTS - 2*obj.N;

if obj.TW 
    MULTS((MIentries_done+1):end,:) = [sabs(isc) obj.dim_chopped*ones(sum(isc),1)];
    VALSUBS(Ventries_done+1:end,:) = [M.subs(is,1) MIentries_done+cumsum(isc)];
    VAL(Ventries_done+1:end) = M.vals(is);
end


[MULTS,mids] = sortrows(sort(MULTS,2));
IC = any(diff([zeros(1,2); MULTS]),2);
MULTS = MULTS(IC,:); 
IC = cumsum(IC);
MULTS = MULTS(:);
N.ind = sparse([1:length(MULTS)/2 1:length(MULTS)/2],MULTS,ones(1,length(MULTS)),length(MULTS)/2,obj.dim_chopped);

midse = zeros(1,length(mids));
midse(mids) = 1:length(mids);
mids = midse;
IC = IC(mids);

VAL(VALSUBS(:,1) == 1) = 0;
VAL(VALSUBS(:,1) == obj.N) = 0;
if strcmp(obj.BC,'sym')
VAL(VALSUBS(:,1) == (1*obj.N+1)) = 0;
VAL(VALSUBS(:,1) == (2*obj.N+1)) = 0;
VAL(VALSUBS(:,1) == (3*obj.N+1)) = 0;
end
for j = 0:2*obj.Q-3 
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+1)) = 0;
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+obj.N)) = 0;
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+obj.N+1)) = 0;
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+2*obj.N)) = 0;
if strcmp(obj.BC,'sym')
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+3*obj.N+1)) = 0;
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+4*obj.N+1)) = 0;
VAL(VALSUBS(:,1) == (4*obj.N+j*bar+5*obj.N+1)) = 0;
end
end
VALSUBS = VALSUBS(VAL ~= 0,:);
VAL = VAL(VAL ~= 0);
sebs = IC(VALSUBS(:,2)); 
[sebs,mic] = sort(sebs);
ice = diff([0; sebs]) ~= 0;
N.ind = N.ind(sebs(ice),:);
ic = cumsum(ice);
n_I = size(N.ind,1);
N.coeffs = sparse(VALSUBS(mic,1),ic,-VAL(mic),obj.dim_chopped,n_I);
% disp(size(N.coeffs))
% disp(size(N.ind)) 

end

function [multi_inds] = multi_ind_parser(mults_init,startpos,entries_done)
multi_inds = ones(size(mults_init,1),1)*startpos + mults_init;
end

function [vels,vels_subs,multi_inds] = vels_parser(dim,vals,valst,stpos,MIentries_done,vals_subs,mults_init,varargin)

if ~ischar(varargin{1})
    imre = varargin{3};
    con1 = [1 1 -1 -1];
    con2 = [1 -1 1 -1];
    vals1 = vals;
    vals1t = valst;
    vals2 = varargin{1};
    vals2t = varargin{2};
else
    imre = varargin{1};
    spec = varargin{2};
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
    case 'double2'
        con1 = [1 1 -1 -1];
        con2 = [1 -1 1 -1];
end
end
numl = length(vals);
numle = size(mults_init,1);
vels = zeros(4*numl,1);
multi_inds = zeros(4*numle,2);
vels_subs = zeros(4*numl,2);
switch imre
    case 'im'
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = -con(2)*vals;  
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = -con(3)*valst; 
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con(4)*vals; 
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
    case 're' 
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = con(2)*vals;
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = con(3)*valst;
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con(4)*vals; 
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
    case 'double'
vels(1:numl) = con1(1)*vals1+con2(1)*vals2t;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = -con1(2)*vals1-con2(3)*vals2t;
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = -con1(3)*vals1t-con2(2)*vals2;
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con1(4)*vals1-con2(4)*vals2t;
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
    case 'double2'
vels(1:numl) = con1(1)*vals+con2(1)*valst; 
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = con1(2)*vals+con2(3)*valst;
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = con1(3)*valst+con2(2)*vals; 
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con1(4)*vals-con2(4)*valst;
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
end
end

function [vels,vels_subs,multi_inds] = vels_parser_n(dim,vals,valst,stpos,MIentries_done,vals_subs,mults_init,varargin)
if ~ischar(varargin{1})
    imre = varargin{3};
    con1 = [1 1 -1 -1];
    con2 = [1 -1 1 -1];
    vals1 = vals;
    vals1t = valst;
    vals2 = varargin{1};
    vals2t = varargin{2};
else
    imre = varargin{1};
    spec = varargin{2};
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
    case 'double2'
        con1 = [1 1 -1 -1];
        con2 = [1 -1 1 -1];
end
end
numl = length(vals);
numle = size(mults_init,1);
vels = zeros(2*numl,1);
multi_inds = zeros(2*numle,2);
vels_subs = zeros(2*numl,2);
% stpos(2) = stpos(2) + 2*obj.N;
switch imre
    case 'im'
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:end) = -con(2)*vals;  
multi_inds(numle+1:end,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:end,:) = vals_subs + [stpos(1) MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(2*numl+1:3*numl) = -con(3)*zeros(size(valst)); 
% multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
% vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(3*numl+1:end) = -con(4)*zeros(size(vals)); 
% multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
% vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
    case 're'
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:end) = con(2)*vals;
multi_inds(numle+1:end,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:end,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(2*numl+1:3*numl) = con(3)*zeros(size(valst));
% multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
% vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(3*numl+1:end) = -con(4)*zeros(size(vals)); 
% multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
% vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1) MIentries_done];
end
end

function [vels,vels_subs,multi_inds] = vels_parser_q(dim,vals,valst,stpos,MIentries_done,vals_subs,mults_init,varargin)
if ~ischar(varargin{1})
    imre = varargin{3};
    con1 = [1 1 -1 -1];
    con2 = [1 -1 1 -1];
    vals1 = vals;
    vals1t = valst;
    vals2 = varargin{1};
    vals2t = varargin{2};
else
    imre = varargin{1};
    spec = varargin{2};
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
    case 'double2'
        con1 = [1 1 -1 -1];
        con2 = [1 -1 1 -1];
end
end
numl = length(vals);
numle = size(mults_init,1);
vels = zeros(4*numl,1);
multi_inds = zeros(4*numle,2);
vels_subs = zeros(4*numl,2);
% stpos(2) = stpos(2) + 2*obj.N;
switch imre
    case 'im'
vels(1:numl) = con(1)*zeros(size(vals));
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = -con(2)*vals;  
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = -con(3)*valst; 
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con(4)*zeros(size(vals)); 
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
    case 're'
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numl+1:2*numl) = con(2)*zeros(size(vals));
multi_inds(numle+1:2*numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:2*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(2*numl+1:3*numl) = con(3)*zeros(size(valst));
multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(3*numl+1:end) = -con(4)*vals; 
multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1) MIentries_done];
end
end

function [vels,vels_subs,multi_inds] = vels_parser_uniq(dim,vals,valst,stpos,MIentries_done,vals_subs,mults_init,varargin)

if ~ischar(varargin{1})
    imre = varargin{3};
    con1 = [1 1 -1 -1];
    con2 = [1 -1 1 -1];
    vals1 = vals;
    vals1t = valst;
    vals2 = varargin{1};
    vals2t = varargin{2};
else
    imre = varargin{1};
    spec = varargin{2};
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
end
end
valsh = varargin{3};
vals_subsh = varargin{4};
multsh = varargin{5};
numl = length(vals);
numle = size(mults_init,1);
numla = size(multsh,1);
numlh = length(valsh);
vels = zeros(2*numl+2*numlh,1);
multi_inds = zeros(2*numle+2*numla,2);
vels_subs = zeros(2*numl+2*numlh,2);
switch imre
    case 'im'
vels(1:numlh) = con(1)*valsh;
multi_inds(1:numla,:) = multi_ind_parser(multsh,[stpos(2:3)],MIentries_done);
vels_subs(1:numlh,:) = vals_subsh + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numla;
vels(numlh+1:numlh+numl) = -con(2)*vals;  
multi_inds(numla+1:numla+numle,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numlh+1:numlh+numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numlh+numl+1:numlh+2*numl) = -con(3)*valst; 
multi_inds(numla+numle+1:numla+2*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(numlh+numl+1:numlh+2*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numlh+2*numl+1:end) = -con(4)*valsh; 
multi_inds(numla+2*numle+1:end,:) = multi_ind_parser(multsh,[stpos(2:3)+dim],MIentries_done);
vels_subs(numlh+2*numl+1:end,:) = vals_subsh + [stpos(1)+dim MIentries_done]; 
    case 're' 
vels(1:numlh) = con(1)*valsh;
multi_inds(1:numla,:) = multi_ind_parser(multsh,[stpos(2:3)],MIentries_done);
vels_subs(1:numlh,:) = vals_subsh + [stpos(1) MIentries_done]; 
MIentries_done = MIentries_done + numla;
vels(numlh+1:numlh+numl) = con(2)*vals;
multi_inds(numla+1:numle+numla,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numlh+1:numlh+numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numlh+numl+1:numlh+2*numl) = con(3)*valst;
multi_inds(numle+numla+1:2*numle+numla,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
vels_subs(numlh+numl+1:numlh+2*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
MIentries_done = MIentries_done + numle;
vels(numlh+2*numl+1:end) = -con(4)*valsh; 
multi_inds(2*numle+numla+1:end,:) = multi_ind_parser(multsh,[stpos(2:3)+dim],MIentries_done);
vels_subs(numlh+2*numl+1:end,:) = vals_subsh + [stpos(1) MIentries_done]; 
end
end

function [vels,vels_subs,multi_inds] = vels_parser_nq(dim,vals,valst,stpos,MIentries_done,vals_subs,mults_init,varargin)
if ~ischar(varargin{1})
    imre = varargin{3};
    con1 = [1 1 -1 -1];
    con2 = [1 -1 1 -1];
    vals1 = vals;
    vals1t = valst;
    vals2 = varargin{1};
    vals2t = varargin{2};
else
    imre = varargin{1};
    spec = varargin{2};
switch spec
    case 'conj1'
        con = [1 1 -1 -1];
    case 'conj2'
        con = [1 -1 1 -1];
    case 'std'
        con = ones(1,4);
    case 'double2'
        con1 = [1 1 -1 -1];
        con2 = [1 -1 1 -1];
end
end
numl = length(vals);
numle = size(mults_init,1);
vels = zeros(1*numl,1);
multi_inds = zeros(1*numle,2);
vels_subs = zeros(1*numl,2);
% stpos(2) = stpos(2) + 2*obj.N;
switch imre
    case 'im'
% vels(1:numl) = con(1)*vals;
% multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
% vels_subs(1:numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
% MIentries_done = MIentries_done + numle;
vels(numl+1:end) = -con(2)*vals;  
multi_inds(numle+1:end,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
vels_subs(numl+1:end,:) = vals_subs + [stpos(1) MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(2*numl+1:3*numl) = -con(3)*zeros(size(valst)); 
% multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
% vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1) MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(3*numl+1:end) = -con(4)*zeros(size(vals)); 
% multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
% vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
    case 're'
vels(1:numl) = con(1)*vals;
multi_inds(1:numle,:) = multi_ind_parser(mults_init,[stpos(2:3)],MIentries_done);
vels_subs(1:numl,:) = vals_subs + [stpos(1) MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(numl+1:end) = con(2)*vals;
% multi_inds(numle+1:end,:) = multi_ind_parser(mults_init,[stpos(2) stpos(3)+dim],MIentries_done);
% vels_subs(numl+1:end,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(2*numl+1:3*numl) = con(3)*zeros(size(valst));
% multi_inds(2*numle+1:3*numle,:) = multi_ind_parser(mults_init,[stpos(3) stpos(2)+dim],MIentries_done);
% vels_subs(2*numl+1:3*numl,:) = vals_subs + [stpos(1)+dim MIentries_done]; 
% MIentries_done = MIentries_done + numle;
% vels(3*numl+1:end) = -con(4)*zeros(size(vals)); 
% multi_inds(3*numle+1:4*numle,:) = multi_ind_parser(mults_init,[stpos(2:3)+dim],MIentries_done);
% vels_subs(3*numl+1:4*numl,:) = vals_subs + [stpos(1) MIentries_done];
end
end