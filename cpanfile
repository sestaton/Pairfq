requires 'DBD::SQLite',         '1.44';
requires 'Tie::Hash::DBD',      '0.13';
requires 'autodie';

on 'test' => sub {
   requires 'Test::More',       '0.96';
};
