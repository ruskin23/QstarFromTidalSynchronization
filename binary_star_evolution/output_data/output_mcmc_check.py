nohup: ignoring input
2019-02-07 10:44:23,377 INFO sqlalchemy.engine.base.Engine SELECT CAST('test plain returns' AS VARCHAR(60)) AS anon_1
2019-02-07 10:44:23,377 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,377 INFO sqlalchemy.engine.base.Engine SELECT CAST('test unicode returns' AS VARCHAR(60)) AS anon_1
2019-02-07 10:44:23,377 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,378 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("interpolators")
2019-02-07 10:44:23,378 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("model_suites")
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("quantities")
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_dependent_variables")
2019-02-07 10:44:23,379 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("interpolation_parameters")
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("tracks")
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_grids")
2019-02-07 10:44:23,380 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,381 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("interpolator_tracks")
2019-02-07 10:44:23,381 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,381 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_age_nodes")
2019-02-07 10:44:23,381 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_feh_nodes")
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_mass_nodes")
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine PRAGMA table_info("varchange_dependent_values")
2019-02-07 10:44:23,382 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,390 INFO sqlalchemy.engine.base.Engine BEGIN (implicit)
2019-02-07 10:44:23,390 INFO sqlalchemy.engine.base.Engine SELECT quantities.id AS quantities_id, quantities.name AS quantities_name 
FROM quantities
2019-02-07 10:44:23,390 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,391 INFO sqlalchemy.engine.base.Engine SELECT model_suites.id AS model_suites_id, model_suites.name AS model_suites_name 
FROM model_suites
2019-02-07 10:44:23,391 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,393 INFO sqlalchemy.engine.base.Engine SELECT count(*) AS count_1 
FROM (SELECT tracks.id AS tracks_id, tracks.filename AS tracks_filename, tracks.mass AS tracks_mass, tracks.feh AS tracks_feh, tracks.model_suite_id AS tracks_model_suite_id, tracks.checksum AS tracks_checksum 
FROM tracks) AS anon_1
2019-02-07 10:44:23,393 INFO sqlalchemy.engine.base.Engine ()
2019-02-07 10:44:23,393 INFO sqlalchemy.engine.base.Engine COMMIT
2019-02-07 10:44:23,394 INFO sqlalchemy.engine.base.Engine BEGIN (implicit)
2019-02-07 10:44:23,395 INFO sqlalchemy.engine.base.Engine SELECT interpolators.id AS interpolators_id, interpolators.name AS interpolators_name, interpolators.filename AS interpolators_filename, interpolators.checksum AS interpolators_checksum 
FROM interpolators 
WHERE interpolators.name = ?
2019-02-07 10:44:23,395 INFO sqlalchemy.engine.base.Engine ('default',)
2019-02-07 10:44:23,500 INFO sqlalchemy.engine.base.Engine SELECT tracks.id AS tracks_id, tracks.filename AS tracks_filename, tracks.mass AS tracks_mass, tracks.feh AS tracks_feh, tracks.model_suite_id AS tracks_model_suite_id, tracks.checksum AS tracks_checksum 
FROM tracks, interpolator_tracks 
WHERE ? = interpolator_tracks.interpolator_id AND tracks.id = interpolator_tracks.track_id
2019-02-07 10:44:23,501 INFO sqlalchemy.engine.base.Engine (1,)
/usr/local/lib/python3.6/dist-packages/sqlalchemy/sql/sqltypes.py:603: SAWarning: Dialect sqlite+pysqlite does *not* support Decimal objects natively, and SQLAlchemy must convert from floating point - rounding errors and other issues may occur. Please consider storing Decimal numbers as strings or integers on this platform for lossless storage.
  'storage.' % (dialect.name, dialect.driver))
2019-02-07 10:44:23,503 INFO sqlalchemy.engine.base.Engine SELECT model_suites.id AS model_suites_id, model_suites.name AS model_suites_name 
FROM model_suites 
WHERE model_suites.id = ?
2019-02-07 10:44:23,503 INFO sqlalchemy.engine.base.Engine (1,)
2019-02-07 10:44:23,504 INFO sqlalchemy.engine.base.Engine SELECT interpolation_parameters.interpolator_id AS interpolation_parameters_interpolator_id, interpolation_parameters.quantity_id AS interpolation_parameters_quantity_id, interpolation_parameters.nodes AS interpolation_parameters_nodes, interpolation_parameters.smoothing AS interpolation_parameters_smoothing, interpolation_parameters.vs_log_age AS interpolation_parameters_vs_log_age, interpolation_parameters.log_quantity AS interpolation_parameters_log_quantity 
FROM interpolation_parameters 
WHERE ? = interpolation_parameters.interpolator_id
2019-02-07 10:44:23,504 INFO sqlalchemy.engine.base.Engine (1,)
2019-02-07 10:44:26,298 INFO sqlalchemy.engine.base.Engine SELECT varchange_grids.id AS varchange_grids_id 
FROM varchange_grids 
WHERE varchange_grids.name = ? AND varchange_grids.interpolator_id = ?
2019-02-07 10:44:26,298 INFO sqlalchemy.engine.base.Engine ('default', 1)
2019-02-07 10:44:26,299 INFO sqlalchemy.engine.base.Engine SELECT varchange_feh_nodes.value AS varchange_feh_nodes_value 
FROM varchange_feh_nodes 
WHERE varchange_feh_nodes.grid_id = ? ORDER BY varchange_feh_nodes."index"
2019-02-07 10:44:26,299 INFO sqlalchemy.engine.base.Engine (1,)
2019-02-07 10:44:26,299 INFO sqlalchemy.engine.base.Engine SELECT varchange_mass_nodes.value AS varchange_mass_nodes_value 
FROM varchange_mass_nodes 
WHERE varchange_mass_nodes.grid_id = ? ORDER BY varchange_mass_nodes."index"
2019-02-07 10:44:26,299 INFO sqlalchemy.engine.base.Engine (1,)
2019-02-07 10:44:26,300 INFO sqlalchemy.engine.base.Engine SELECT varchange_age_nodes.value AS varchange_age_nodes_value 
FROM varchange_age_nodes 
WHERE varchange_age_nodes.grid_id = ? ORDER BY varchange_age_nodes."index"
2019-02-07 10:44:26,300 INFO sqlalchemy.engine.base.Engine (1,)
2019-02-07 10:44:26,302 INFO sqlalchemy.engine.base.Engine COMMIT
calculatingforthevalues = 4.207876469222957e-08	6.850758421971131

FINISHED PLANET-STAR EVOLUTION
send_angmom =  [1.1295307  0.00702629]
Target parameters: 
age =  1.6835201199819239
convective phase lag =  4.207876469222957e-08
self.disk_lock_frequency  =  4.570685778068612
Trying P0 = 5.2663825, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.323207351536312, P* = 5.323207351536315
Trying P0 = 2.63319125, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 2.1897340385596147, P* = 2.1897340385596156
For Pdisk = 1.37467015066492, orbital period range: 2.63319125 < Porb < 5.2663825
Trying P0 = 2.63319125, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 2.1904416683566312, P* = 2.1904416683566317
Trying P0 = 5.2663825, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.333796721053445, P* = 5.333796721053447
Trying P0 = 5.209909548638301, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.267682666223794, P* = 5.267682666223796
Trying P0 = 5.208799448802316, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.266591095663324, P* = 5.266591095663326
Trying P0 = 5.208587328090274, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.266382515448898, P* = 5.2663825154489
Trying P0 = 5.20858422379661, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.2663794649209095, P* = 5.266379464920912
Trying P0 = 5.208587328090274, Pdisk = 1.37467015066492
BINARY CONFIGURATION COMPLETE
BINARY EVOLUTION COMPLETE
Got Porb = 5.266382515448898, P* = 5.2663825154489
183.92user 2.00system 3:03.79elapsed 101%CPU (0avgtext+0avgdata 209068maxresident)k
0inputs+56outputs (0major+108508minor)pagefaults 0swaps
