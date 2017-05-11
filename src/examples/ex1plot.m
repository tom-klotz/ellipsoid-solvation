A = importdata('out/ex1.txt');
x = 1:18;
figure;
semilogy(x,A(:,3));
hold on;
semilogy(x,A(:,4));
xlabel('Max expansion order (N)')
ylabel('Absolute Error')
legend('Ellipsoidal', 'Spherical')