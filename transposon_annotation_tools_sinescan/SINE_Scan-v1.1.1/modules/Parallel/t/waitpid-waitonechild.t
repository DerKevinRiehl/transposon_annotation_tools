use strict;
use warnings;

use Test::More;
use Test::Warn;

use Parallel::ForkManager;

my $pm = Parallel::ForkManager->new(4);

local $SIG{ALRM} = sub {
    fail "test hanging, forever waiting for child process";
    exit 1;
};

for ( 1 ) {
    $pm->start and last;
    sleep 2;
    $pm->finish;
}

my $pid = waitpid -1, 0;

diag "code outside of P::FM stole $pid";

TODO: {
    local $TODO = 'MacOS and FreeBDS seem to have issues with this';

    eval {
        alarm 10;
        warning_like {
            $pm->wait_one_child;
        } qr/child process '\d+' disappeared. A call to `waitpid` outside of Parallel::ForkManager might have reaped it\./,
            "got the missing child warning";
        pass "wait_one_child terminated";
    };

    is $pm->running_procs => 0, "all children are accounted for";

}

done_testing;
