% Completed by Matthew Goldberg on Apr 30 2018
% Unit tests for math thesis functions
% How to use: 
% testCase1 = ThesisTest;
% result1 = run(testCase1);
% disp(result1)
classdef ThesisTest < matlab.unittest.TestCase
    % Basic unit test for functions used in math thesis
    
    properties
    end
    
    methods (Test)
        function testTridag(testCase)
            actualx = tridiag_solve([0 2 1 1],[6 4 4 6],[1 1 2 0],[8 13 22 27]);
            expectedx = [1 2 3 4];
            testCase.verifyEqual(actualx,expectedx,'AbsTol',1e-6);
        end

        function testFEM(testCase)
            Uanalytical = @(x,t)(exp(-t)).*sin(pi*x) + (exp(-4*t)).*sin(2*pi*x);
            MyU = Uanalytical(linspace(0,1,33)',linspace(0,1,33));
            [actualU,~] = MyDiffusionFEM(33,33,pi^-2,1,1,@(x) sin(pi*x)+sin(2*pi*x),Uanalytical);
            
            % Method 1:
            % expectedU = actualU; % this is silly, but in practice we
            % would have some analytical solution in mind -> see MyU
            
            % Method 2:
            % Are the mean values of these U's within 1e-3 of one another?
            % This is actually a measure of the error in our calculation!
            testCase.verifyEqual(mean(actualU(:)),mean(MyU(:)),'AbsTol',1e-3);
        end
    end
    
end