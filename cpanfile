requires 'BerkeleyDB',          '>= 0.54';
requires 'IPC::System::Simple', '>= 1.21';
requires 'List::MoreUtils',     '>= 0.33';

on 'test' => sub {
   requires 'Test::More', '>= 0.96';
};
