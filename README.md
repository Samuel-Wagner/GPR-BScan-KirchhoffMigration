# GPR-BScan-KirchhoffMigration
This is a script that shows how to migrate a ground penetrating radar (GPR) B-scan. Migration is one of several imaging algorithms. 
This is an implementation of 2D Kirchhoff Migration. The scans have been mostly pre-processed already. To run, run main.m. 

The real important scripts are:
1. utils/Bscan_migration_v3.m (actual Kirchhoff Migration implementation)
2. utils/GPR_transmission_angles_v4.m (quick-solver for GPR transmission angles)

The rest are useful, nice scripts that help process data or help with display.

To run, run main.m and read the comments there to change the data source/etc.
There is an example of migration using a "Scan" object (what I personally like to use) and migration using manual data. I also include a timing test - this is a faster migration algorithm than those I've previously shown online.

Sample of images generated with main.m:

Raw and background-removed data (BKGR)
![Raw and background-removed data](https://github.com/Samuel-Wagner/GPR-BScan-KirchhoffMigration/blob/main/imgs/fig1.png)
Migration before and after a matched-filter
![Migrated normal and matched-filtered data](https://github.com/Samuel-Wagner/GPR-BScan-KirchhoffMigration/blob/main/imgs/fig3.png)
Horizontal and Vertical Resolution Slices of Target 2 with and without matched-filter
![Horizontal and Vertical Target Slices](https://github.com/Samuel-Wagner/GPR-BScan-KirchhoffMigration/blob/main/imgs/fig4.png)
Timing test versus pixel width (40cm x 1.3m calculation domain)
![Timing test](https://github.com/Samuel-Wagner/GPR-BScan-KirchhoffMigration/blob/main/imgs/fig5.png)
