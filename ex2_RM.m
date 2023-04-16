%% Hopfield Model

P=15; %number of images
N=100; %size of cells in image matrix
K=50; %max percent of noise added to the image
it=3500; %iterations
E=30; %number of repititions of net tests for the mean performance

perf=zeros(P,K/5); %create the matrix for a plot later. jumps of 5 in noise percentage

clear Image
%generate small p random images
for m=1:P
    Image{m}=randi([0, 1], [sqrt(N), sqrt(N)]);
    for i=1:N
        if Image{m}(i)==0
            Image{m}(i)=-1;
        end
    end
end

for p=1:P   
    %build the net connections
    W=zeros(N,N); %net's weights
    for i=1:N
        for j=1:N
            if i~=j
                for m=1:p
                    Wtemp=(1/N)*Image{m}(i)*Image{m}(j);
                    W(i,j)=W(i,j)+Wtemp;
                end
            else
                W(i,j)=0;
            end
        end
    end
    
    %add noise to images
    counter=0;
    for k=5:5:K
        counter=counter+1;
        for m=1:p
            noisecells=randsample(N,N*k/100);
            Image_Noise{m}=Image{m};
            for i=1:length(noisecells)
                Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
            end
        end
        
        for e=1:E %number of repitions for mean performance
            Image_Pick=randsample(p,1);
            Net=Image_Noise{Image_Pick}; %choose Image to present to the net
            Cor=0;
            t=1;
            conv(e)=0;
            %update the net one neuron at a a time
            for t=1:it
                up=randsample(N,1); %random update index
                Net(up)=0;
                for i=1:N
                    Net(up)=Net(up)+Net(i)*W(up,i);
                end
                Net(up)=sign(Net(up));
                r=corrcoef(Net,Image{Image_Pick}); %find the correlation between original image and the image with noise
                Cor(t)=r(2);
                if Cor(t)==1 && conv(e)==0
                    conv(e)=t; %iteration of convergence
                end
                if t>100 && Cor(t)==1 && Cor(t-100)==1
                    break
                end
            end
            
            if Cor(t)==1
                perf(p,counter)=perf(p,counter)+1/E;
            else
                conv(e)=it;
            end
        end
        mean_conv(p,counter)=mean(conv(conv~=0)); %number of steps for convergence on average
    end
end

%plot the results of the different nets
figure();
for i=1:P
    plot(perf(i,:),'color',[(i/P) 0 1-(i/P)],'LineWidth',2);
    hold on
end
xlabel('Noise Percentage');
ylabel('Successful Convergence Percentage');
title('Performances of Nets: Succession rates');
legend(['Blue is 1 memory, Red is ',num2str(P),' Memories'],'Location','southwest');
ylim([0 1]);
yticklabels({'0%' '10%' '20%' '30%' '40%' '50%' '60%' '70%' '80%' '90%' '100%'});
xlim([1 K/5]);
xticklabels({'5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%' '45%' '50%'});

figure();
a=perf>=0.99;
b=mean_conv;
b(~a)=NaN;
for i=1:P
    plot(b(i,:),'color',[0 (i/P) 0],'LineWidth',2);
    hold on
end
xlabel('Noise Percentage');
ylabel('Average Convergence Iterations');
title('Performances of Nets: Convergence Speed');
legend(['dark green is 1 memory, light green is ',num2str(P),' Memories'],'Location','northwest');
% ylim([0 it]);
xlim([1 K/5]);
xlab=['5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%'];
xticklabels({'5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%' '45%' '50%'});

%% single net testing:
P=4; %number of images
N=100; %size of cells in image matrix
K=20; %max percent of noise added to the image
it=2500; %iterations

%random images
for m=1:P
    Image{m}=randi([0, 1], [sqrt(N), sqrt(N)]);
    for i=1:N
        if Image{m}(i)==0
            Image{m}(i)=-1; 
        end
    end
end
W=zeros(N,N); %net's weights
for i=1:N
    for j=1:N
        if i~=j
            for m=1:P
                Wtemp=(1/N)*Image{m}(i)*Image{m}(j);
                W(i,j)=W(i,j)+Wtemp;
            end
        else
            W(i,j)=0;
        end
    end
end

% % check that the net is stable with one memory
% Net=Image{1}; %choose Image to present to the net
% Cor=0;
% t=1;
% for t=1:it
%     up=randsample(N,1); %update index
%     Net(up)=0;
%     for i=1:N
%         Net(up)=Net(up)+Net(i)*W(up,i);
%     end
%     Net(up)=sign(Net(up));
%     r=corrcoef(Net,Image{1});
%     Cor(t)=r(2);
% end
% figure()
% plot(Cor);
% xlabel('Iterations');
% ylabel('Correlation (R^2)');
% title('Stability of the Net');
% legend(['If y=1, net is stable']);

% add noise to images and present to the net
for m=1:P
    noisecells=randsample(N,N*K/100);
    Image_Noise{m}=Image{m};
    for i=1:length(noisecells)
        Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
    end
end
Net=Image_Noise{1}; %choose Image to present to the net
Cor=0;
t=1;
conv=0;
for t=1:it
    up=randsample(N,1); %update index
    Net(up)=0;
    for i=1:N
        Net(up)=Net(up)+Net(i)*W(up,i);
    end
    Net(up)=sign(Net(up));
    r=corrcoef(Net,Image{1});
    Cor(t)=r(2);
    if Cor(t)==1 && conv==0
        conv=t; %iteration of convergence
    end
    if t>100 && Cor(t)==1 && Cor(t-100)==1
        break
    end
end
if t<it
    disp(['stopped updating 100 steps after convergence (', num2str(t), ' iterations total)' ])
else
    disp(['Net did not reach target memory after ', num2str(it), ' iterations'])
end
figure()
plot(Cor);
xlabel('Iterations');
ylabel('Correlation (R^2)');
title('Net progress');
if Cor(end)==1
    legend(['Convergence in ',num2str(conv),' Iterations']);
else
    legend(['Net did not reach target memory']);
end
ylim([0 1]);
xlim([0 it]);
%% Q1 & Q3

%for question 1 we will pick 5 memories and test them with different noise
%levels. 

P=5; %number of images
N=100; %size of cells in image matrix
K=50; %max percent of noise added to the image
it=3500; %iterations
E=30; %number of repititions of net tests for the mean performance

p=P;
clear Image
clear mean_conv
perf=zeros(1,K/5);
%generate small p random images
for m=1:p
    Image{m}=randi([0, 1], [sqrt(N), sqrt(N)]);
    for i=1:N
        if Image{m}(i)==0
            Image{m}(i)=-1;
        end
    end
end

%build the net connections
W=zeros(N,N); %net's weights
for i=1:N
    for j=1:N
        if i~=j
            for m=1:p
                Wtemp=(1/N)*Image{m}(i)*Image{m}(j);
                W(i,j)=W(i,j)+Wtemp;
            end
        else
            W(i,j)=0;
        end
    end
end

%add noise to images
counter=0;
for k=5:5:K
    counter=counter+1;
    for m=1:p
        noisecells=randsample(N,N*k/100);
        Image_Noise{m}=Image{m};
        for i=1:length(noisecells)
            Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
        end
    end
    
    for e=1:E %number of repitions for mean performance
        Image_Pick=randsample(p,1);
        Net=Image_Noise{Image_Pick}; %choose Image to present to the net
        Cor=0;
        t=1;
        conv(e)=0;
        %update the net one neuron at a a time
        for t=1:it
            up=randsample(N,1); %random update index
            Net(up)=0;
            for i=1:N
                Net(up)=Net(up)+Net(i)*W(up,i);
            end
            Net(up)=sign(Net(up));
            r=corrcoef(Net,Image{Image_Pick}); %find the correlation between original image and the net
            Cor(t)=r(2);
            if Cor(t)==1 && conv(e)==0
                conv(e)=t; %iteration of convergence
            end
            if t>100 && Cor(t)==1 && Cor(t-100)==1
                break
            end
        end
        
        if Cor(t)==1
            perf(1,counter)=perf(1,counter)+1/E;
        else
            conv(e)=it;
        end
    end
    mean_conv(1,counter)=mean(conv(conv~=0)); %number of steps for convergence on average
end


%plot the results of the different nets
figure();
plot(perf,'LineWidth',2);
hold on
xlabel('Noise Percentage');
ylabel('Successful Convergence Percentage');
title([num2str(P) ' Memories: Succesful Convergence vs. Noise']);
ylim([0 1]);
yticklabels({'0%' '10%' '20%' '30%' '40%' '50%' '60%' '70%' '80%' '90%' '100%'});
xlim([1 K/5]);
xticklabels({'5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%' '45%' '50%'});

figure();
plot(mean_conv,'LineWidth',2);
hold on
xlabel('Noise Percentage');
ylabel('Average Convergence Iterations');
title([num2str(P) ' Memories: Iteration vs. Noise']);
% ylim([0 it]);
a=find(perf>=0.99);
xticklabels({'5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%' '45%' '50%'});
if length(a)>=1
    xlim([1 a(end)]);
end

%% Q2

%for Q2 we will choos 1-20 memories, 15% noise

P=20; %number of images
N=100; %size of cells in image matrix
K=30; %max percent of noise added to the image
it=3500; %iterations
E=30; %number of repititions of net tests for the mean performance

perf=zeros(P,1); %create the matrix for a plot later. jumps of 5 in noise percentage

clear Image
%generate small p random images
for m=1:P
    Image{m}=randi([0, 1], [sqrt(N), sqrt(N)]);
    for i=1:N
        if Image{m}(i)==0
            Image{m}(i)=-1;
        end
    end
end

for p=1:P
    %build the net connections
    W=zeros(N,N); %net's weights
    for i=1:N
        for j=1:N
            if i~=j
                for m=1:p
                    Wtemp=(1/N)*Image{m}(i)*Image{m}(j);
                    W(i,j)=W(i,j)+Wtemp;
                end
            else
                W(i,j)=0;
            end
        end
    end
    
    %add noise to images
    for m=1:p
        noisecells=randsample(N,N*K/100);
        Image_Noise{m}=Image{m};
        for i=1:length(noisecells)
            Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
        end
    end
    
    for e=1:E %number of repitions for mean performance
        Image_Pick=randsample(p,1);
        Net=Image_Noise{Image_Pick}; %choose Image to present to the net
        Cor=0;
        t=1;
        conv(e)=0;
        %update the net one neuron at a a time
        for t=1:it
            up=randsample(N,1); %random update index
            Net(up)=0;
            for i=1:N
                Net(up)=Net(up)+Net(i)*W(up,i);
            end
            Net(up)=sign(Net(up));
            r=corrcoef(Net,Image{Image_Pick}); %find the correlation between original image and the net
            Cor(t)=r(2);
            if Cor(t)==1 && conv(e)==0
                conv(e)=t; %iteration of convergence
            end
            if t>100 && Cor(t)==1 && Cor(t-100)==1
                break
            end
        end
        
        if Cor(t)==1
            perf(p)=perf(p)+1/E;
        else
            conv(e)=it;
        end
    end
    mean_conv(p)=mean(conv(conv~=0)); %number of steps for convergence on average
end

%plot the results of the net
figure();
plot(perf,'LineWidth',2);
hold on
xlabel('Number of Memories');
ylabel('Successful Convergence Percentage');
title([num2str(K) '% Noise: Memories vs. Successful Convergence']);
ylim([0 1]);
yticklabels({'0%' '10%' '20%' '30%' '40%' '50%' '60%' '70%' '80%' '90%' '100%'});
xlim([1 20]);

%% Q4

%we have seen convergence for random noise, let's check deterministic noise
clear Image
p=7; %number of images
N=100; %size of cells in image matrix
K=50; %max percent of noise added to the image
it=3500; %iterations
E=30; %number of repititions of net tests for the mean performance
clear mean_conv
%generate small p random images
for m=1:p
    Image{m}=randi([0, 1], [sqrt(N), sqrt(N)]);
    for i=1:N
        if Image{m}(i)==0
            Image{m}(i)=-1;
        end
    end
end
P=p;
clear p
counter2=0;
for p=3:2:P
    counter2=counter2+1;
    %build the net connections
    W=zeros(N,N); %net's weights
    for i=1:N
        for j=1:N
            if i~=j
                for m=1:p
                    Wtemp=(1/N)*Image{m}(i)*Image{m}(j);
                    W(i,j)=W(i,j)+Wtemp;
                end
            else
                W(i,j)=0;
            end
        end
    end
    
    %add noise to images
    counter=0;
    for k=5:5:K
        counter=counter+1;
        for m=1:p
            noisecells=(1:k);
            Image_Noise{m}=Image{m};
            for i=1:length(noisecells)
                Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
            end
        end
        
        for e=1:E %number of repitions for mean performance
            Image_Pick=randsample(p,1);
            Net=Image_Noise{Image_Pick}; %choose Image to present to the net
            Cor=0;
            t=1;
            conv(e)=0;
            %update the net one neuron at a a time
            for t=1:it
                up=randsample(N,1); %random update index
                Net(up)=0;
                for i=1:N
                    Net(up)=Net(up)+Net(i)*W(up,i);
                end
                Net(up)=sign(Net(up));
                r=corrcoef(Net,Image{Image_Pick}); %find the correlation between original image and the net
                Cor(t)=r(2);
                if Cor(t)==1 && conv(e)==0
                    conv(e)=t; %iteration of convergence
                end
                if t>100 && Cor(t)==1 && Cor(t-100)==1
                    break
                end
            end
        end
        mean_conv_det(counter2,counter)=mean(conv(conv~=0)); %number of steps for convergence on average in determenistic noise
    end
    
    %add random noise to images
    
    counter=0;
    for k=5:5:K
        counter=counter+1;
        for m=1:p
            noisecells=randsample(N,N*k/100);
            Image_Noise{m}=Image{m};
            for i=1:length(noisecells)
                Image_Noise{m}(noisecells(i))=Image_Noise{m}(noisecells(i))*(-1);
            end
        end
        
        for e=1:E %number of repitions for mean performance
            Image_Pick=randsample(p,1);
            Net=Image_Noise{Image_Pick}; %choose Image to present to the net
            Cor=0;
            t=1;
            conv(e)=0;
            %update the net one neuron at a a time
            for t=1:it
                up=randsample(N,1); %random update index
                Net(up)=0;
                for i=1:N
                    Net(up)=Net(up)+Net(i)*W(up,i);
                end
                Net(up)=sign(Net(up));
                r=corrcoef(Net,Image{Image_Pick}); %find the correlation between original image and the net
                Cor(t)=r(2);
                if Cor(t)==1 && conv(e)==0
                    conv(e)=t; %iteration of convergence
                end
                if t>100 && Cor(t)==1 && Cor(t-100)==1
                    break
                end
            end
        end
        mean_conv_rand(counter2,counter)=mean(conv(conv~=0)); %number of steps for convergence on average in determenistic noise
    end
end

figure();
for i=1:size(mean_conv_rand,1)
    plot(mean_conv_det(i,:),'LineWidth',2,'color',[(i/3) 0 1-(i/3)]);
    hold on
    plot(mean_conv_rand(i,:),'--','LineWidth',2,'color',[(i/3) 0 1-(i/3)]);
    hold on
end
xlabel('Noise Percentage');
ylabel('Average Convergence Iterations');
title(['Determenistic Noise vs. Random Noise']);
xticklabels({'5%' '10%' '15%' '20%' '25%' '30%' '35%' '40%' '45%' '50%'});
legend('3 memories, det noise','3 memories, rand noise','5 memories, det noise','5 memories, rand noise','7 memories, det noise','7 memories, rand noise','Location','southeast')


