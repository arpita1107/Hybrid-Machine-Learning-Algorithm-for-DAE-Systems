function NN = electrode_NN(X)

electrode_global;

x1 = X(1);
x2 = X(2);
w1 = X(3);
w2 = X(4);
w3 = X(5);

X1 = [x1;x2];

% Input weights to the hidden layer

W1 = [w1 w2];

% Input to the hidden layer

h1 = W1*X1 + b;

% Activation function - tanh

out = (exp(h1) - exp(-h1))/(exp(h1) + exp(-h1));

% NN output

op = (w3*out)+b;

NN = (exp(op) - exp(-op))/(exp(op) + exp(-op));

end