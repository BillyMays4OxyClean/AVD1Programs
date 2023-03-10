%% Argument validation testing

Test('ea;')
Test('eaguh','eaea')

function Test(blah,blaher)
    if nargin == 1
        fprintf('yo this dun exist\n')
    elseif nargin == 2
        fprintf('aye you did this right\n')
    end
end