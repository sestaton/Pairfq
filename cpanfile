requires 'DBD::SQLite',         '1.44';
requires 'Tie::Hash::DBD',      '0.13';
requires 'IPC::System::Simple', '1.21';
requires 'List::MoreUtils',     '0.33';

on 'test' => sub {
   requires 'Test::More',       '0.96';
};
