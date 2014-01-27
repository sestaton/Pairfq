requires 'IPC::System::Simple', '>= 1.21';
requires 'BerkeleyDB', '>= 0.54';

on 'test' => sub {
   requires 'Test::More', '>= 0.96';
};
