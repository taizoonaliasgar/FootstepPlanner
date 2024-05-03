%clear
HLflytrot = struct;
%HLflytrot.foot_state = readmatrix('foot_state_mit.txt');
HLflytrot.state = readmatrix('state.txt');
HLflytrot.forces = readmatrix('footforce.txt');
%HLflytrot.mpc_state = readmatrix('mpc_state_mit.txt');
HLflytrot.desired_state = readmatrix('desired_state.txt');
HLflytrot.contact_index = readmatrix('contact_index.txt');
HLflytrot.foot_position = readmatrix('foot_position.txt');
HLflytrot.force_initial = readmatrix('forceinitial.txt');

params = {'Linewidth',2};
paramsd = {'Linewidth',4};
paramst = { 'Interpreter', 'Latex', 'FontSize',25};
paramstick = { 'TickLabelInterpreter', 'Latex', 'FontSize',25};

figure(101)
hold off
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,2),'b-.',params{:})
hold on
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,3),'r-.',params{:})
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,4),'k-.',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,1),'bo ',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,2),'ro ',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,3),'ko ',params{:})
set(gca,paramstick{:})
legend('$x$','$y$','$z$',paramst{:});%'$x_{des}$','$y_{des}$','$z_{des}$',paramst{:})
title('HLMPC CoM coordinates',paramst{:});

figure(102)
hold off
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,8),'b-.',params{:})
hold on
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,9),'r-.',params{:})
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,10),'k-.',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,4),'b.',paramsd{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,5),'r.',paramsd{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,6),'k.',paramsd{:})
set(gca,paramstick{:})
legend('Roll','Pitch','Yaw',paramst{:})
%legend('$\dot{x}$','$\dot{y}$','$\dot{z}$','$\dot{x}_{des}$','$\dot{y}_{des}$','$\dot{z}_{des}$',paramst{:})
title('HLMPC trunk orientation',paramst{:});
%title('HLMPC CoM velocities',paramst{:});


% figure(103)
% hold off
% plot(0.001*HLflytrot.state(:,7),HLflytrot.state(:,7),'b',params{:})
% hold on
% plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.state(:,8),'r',params{:})
% plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.state(:,9),'k',params{:})
% set(gca,paramstick{:})
% legend('Roll','Pitch','Yaw',paramst{:})
% title('HLMPC trunk orientation',paramst{:});

figure(104)
hold off
plot(0.01*HLflytrot.state(:,1),HLflytrot.forces(:,[3,6,9,12]),params{:})
legend('Front Right','Front Left','Rear Right','Rear Left',paramst{:})
set(gca,paramstick{:})
title('Z force',paramst{:})


% figure(105)
% hold off
% plot(0.001*HLflytrot.state(:,1),HLflytrot.forces(:,[1,4,7,10]),params{:})
% legend('Front Right','Front Left','Rear Right','Rear Left',paramst{:})
% set(gca,paramstick{:})
% title('X force',paramst{:})

figure(106)
hold off
plot(0.001*HLflytrot.state(:,1),HLflytrot.forces(:,[2,5,8,11]),params{:})
legend('Front Right','Front Left','Rear Right','Rear Left',paramst{:})
set(gca,paramstick{:})
title('Y force',paramst{:})

figure(107)
hold off
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,5),'b-.',params{:})
hold on
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,6),'r-.',params{:})
plot(0.01*HLflytrot.state(:,1),HLflytrot.state(:,7),'k-.',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,1),'bo ',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,2),'ro ',params{:})
%plot(0.001*HLflytrot.foot_state(:,7),HLflytrot.desired_state(:,3),'ko ',params{:})
set(gca,paramstick{:})
legend('$v_x$','$v_y$','$v_z$',paramst{:});%'$x_{des}$','$y_{des}$','$z_{des}$',paramst{:})
title('HLMPC CoM velocity',paramst{:});















