[IM,map] = imread('SuperMario.png');
imshow(IM,map);
hold on;
set(gca,'Ydir','reverse');

Innertrackx = [];
Innertracky = [];
Outertrackx = [];
Outertracky = [];
Midtrackx = [];
Midtracky = [];

% define control points ---------------------------------------------------
% A, 1
p{1} = [542  744];
q{1} = [635  922];
s{1} = [613 696];
v{1} = [702 877];
r{1} = [824 1031];
d{1} = 0;
w{1} = [822 942];
e{1} = 0;
% B, 2
p{2} = [978 950];
q{2} = [978 550];
s{2} = [895 850];
v{2} = [895 550];
r{2} = [947 453];
d{2} = 1;
w{2} = [873 504];
e{2} = 1;
% C, 3
p{3} = [883 417];
q{3} = [341 144];
s{3} = [843 485];
v{3} = [315 219];
r{3} = [134 116];
d{3} = -0.5;
w{3} = [200 190];
e{3} = 0;
% D, 4
p{4} = [40 286];
q{4} = [40 756];
s{4} = [118 290];
v{4} = [118 692];
r{4} = [126 860];
d{4} = 0;
w{4} = [155 759];
e{4} = 0;
% E, 5
p{5} = [229 823];
q{5} = [433 719];
s{5} = [203 742];
v{5} = [401 643];
r{5} = [497 705];
d{5} = 0;
w{5} = [497 612];
e{5} = 0;

for i = 1:5
    d_1{i} = (q{i}(2)-p{i}(2))/(q{i}(1)-p{i}(1));
    d_2{i} = d_1{i};
end
d_1{2} = 4;
d_2{2} = -2;
d_1{4} = 4;
d_2{4} = -4;

for i = 1:5
    e_1{i} = (v{i}(2)-s{i}(2))/(v{i}(1)-s{i}(1));
    e_2{i} = e_1{i};
end
e_1{2} = 4;
e_2{2} = -4;
e_1{4} = 4;
e_2{4} = -4;

N = 10;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
for i = 1:5
    col = 'b';
    % outer part
    [gamma{i}{1}] = splineinter(q{i}(1),q{i}(2),d_1{i},r{i}(1),r{i}(2),d{i});
    [gamma{i}{2}] = splineinter(r{i}(1),r{i}(2),d{i},p{ mod(i,5)+1}(1),p{mod(i,5)+1}(2), d_2{mod(i,5)+1});
    x0 = q{i}(1);
    x_ = linspace(q{i}(1),r{i}(1),N);
    coef = gamma{i}{1};
    plot(x_,coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3,'LineWidth',2,'color',col)
    x_ = x_(2:end);
    Outertrackx = [Outertrackx x_];
    Outertracky = [Outertracky coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3];
    x0 = r{i}(1);
    x_ = linspace(r{i}(1),p{mod(i,5)+1}(1),N);
    coef = gamma{i}{2};
    plot(x_,coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3,'LineWidth',2,'color',col)
    x_ = x_(2:end);
    Outertrackx = [Outertrackx x_];
    Outertracky = [Outertracky (coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3)];
    % inner part
    [theta{i}{1}] = splineinter(v{i}(1),v{i}(2),e_1{i},w{i}(1),w{i}(2),e{i});
    [theta{i}{2}] = splineinter(w{i}(1),w{i}(2),e{i},s{ mod(i,5)+1}(1),s{mod(i,5)+1}(2), e_2{mod(i,5)+1});
    x0 = v{i}(1);
    x_ = linspace(v{i}(1),w{i}(1),N);
    coef = theta{i}{1};
    plot(x_,coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3,'LineWidth',2,'color',col)
    x_ = x_(2:end);
    Innertrackx = [Innertrackx x_];
    Innertracky = [Innertracky coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3];
    x0 = w{i}(1);
    x_ = linspace(w{i}(1),s{mod(i,5)+1}(1),N);
    coef = theta{i}{2};
    plot(x_,coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3,'LineWidth',2,'color',col)
    x_ = x_(2:end);    
    Innertrackx = [Innertrackx x_];
    Innertracky = [Innertracky coef(1) + coef(2)*(x_-x0) + coef(3)*(x_-x0).^2 + coef(4)*(x_-x0).^3];
end
for i = 1:5
    col = 'b';
    line([p{i}(1) q{i}(1)],[p{i}(2) q{i}(2)],'LineWidth',2,'color',col);
    line([s{i}(1) v{i}(1)],[s{i}(2) v{i}(2)],'LineWidth',2,'color',col);
end
for i = 1:length(Innertrackx)
    Midtrackx(i) = 0.5*(max(Outertrackx(i),Innertrackx(i)) + min(Outertrackx(i),Innertrackx(i)));
    Midtracky(i) = 0.5*(max(Outertracky(i),Innertracky(i)) + min(Outertracky(i),Innertracky(i)));
end
Track = [Outertrackx; Outertracky; Innertrackx; Innertracky; Midtrackx; Midtracky];
Finish = [[896 977];[622 622]];
%plot(Outertrackx,Outertracky,'*r');
%plot(Innertrackx,Innertracky,'*r');
%plot(Midtrackx,Midtracky,'*w');

% finish line
line(Finish(1,:),Finish(2,:),'LineWidth',2,'color',[0.5 0.5 0.5]);
axis([0 1024 0 1100]);

xlabel('x'); ylabel('y')
%text(-30,-30,'(x,y)=(0,0)')
%text(-30,1100+30,'(x,y)=(0,1100)')
%text(1024-50,-30,'(x,y)=(1024,0)')