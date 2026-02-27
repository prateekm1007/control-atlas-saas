

# ══════════════════════════════════════════════════════════════════
# PHASE B1: ENHANCED COMPARISON VISUALIZATION
# ══════════════════════════════════════════════════════════════════

# Add this to your existing comparison tab/section
if hasattr(st.session_state, 'comparison_baseline') and st.session_state.comparison_baseline:
    st.markdown("---")
    st.header("📊 Before/After Improvement Analysis")
    
    # Fetch comparison data (reuse existing comparison_engine logic)
    from comparison_engine import compare_audits
    
    try:
        # Load both audits (you'll need to implement audit retrieval by ID)
        # For now, showing structure assuming audits are in session state
        
        baseline_id = st.session_state.comparison_baseline
        refined_id = st.session_state.comparison_refined
        
        # Placeholder: retrieve audits from storage
        st.info(f"Comparing: {baseline_id} → {refined_id}")
        
        # Mock comparison structure (replace with real data)
        comparison = {
            "verdict_change": {
                "baseline": "INDETERMINATE",
                "refined": "PASS",
                "improved": True
            },
            "coverage_delta": +35.2,
            "det_score_delta": +15,
            "laws": {}
        }
        
        # Display high-level metrics
        col1, col2, col3 = st.columns(3)
        
        with col1:
            verdict_color = "🟢" if comparison["verdict_change"]["improved"] else "🔴"
            st.metric(
                "Verdict",
                comparison["verdict_change"]["refined"],
                delta=f"{verdict_color} {comparison['verdict_change']['baseline']} → {comparison['verdict_change']['refined']}"
            )
        
        with col2:
            coverage_delta = comparison["coverage_delta"]
            st.metric(
                "Coverage",
                "75.2%",  # refined value
                delta=f"{coverage_delta:+.1f}%"
            )
        
        with col3:
            det_delta = comparison["det_score_delta"]
            st.metric(
                "Deterministic Score",
                "85/100",  # refined value
                delta=f"{det_delta:+d} points"
            )
        
        st.markdown("---")
        
        # Law-by-law breakdown
        st.subheader("🔬 Law-Level Improvements")
        
        # Example law improvement (replace with real iteration)
        improvements = [
            {
                "law_id": "LAW-150",
                "title": "Rotamer Audit",
                "baseline_status": "VETO",
                "refined_status": "PASS",
                "baseline_obs": 54.55,
                "refined_obs": 8.23,
                "delta": -46.32,
                "fixed_residues": ["A:6", "A:7", "A:8", "A:9"],
                "still_failing": ["A:10", "A:14"]
            },
            {
                "law_id": "LAW-130",
                "title": "Clashscore",
                "baseline_status": "VETO",
                "refined_status": "PASS",
                "baseline_obs": 26340.0,
                "refined_obs": 12.5,
                "delta": -26327.5,
                "fixed_pairs": [("A:1", "A:3")],
                "still_failing": []
            }
        ]
        
        for imp in improvements:
            with st.expander(f"{'✅' if imp['refined_status'] == 'PASS' else '⚠️'} {imp['law_id']}: {imp['title']}", expanded=True):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown(f"**Before:** {imp['baseline_status']}")
                    st.metric("Observed", f"{imp['baseline_obs']:.2f}")
                
                with col2:
                    st.markdown(f"**After:** {imp['refined_status']}")
                    st.metric("Observed", f"{imp['refined_obs']:.2f}", delta=f"{imp['delta']:.2f}")
                
                with col3:
                    improvement_pct = (abs(imp['delta']) / imp['baseline_obs'] * 100) if imp['baseline_obs'] != 0 else 0
                    st.metric("Improvement", f"{improvement_pct:.1f}%")
                
                # Residue-level improvements
                if "fixed_residues" in imp and imp["fixed_residues"]:
                    st.success(f"✅ Fixed residues: {', '.join(imp['fixed_residues'][:10])}")
                
                if "still_failing" in imp and imp["still_failing"]:
                    st.warning(f"⚠️ Still failing: {', '.join(imp['still_failing'][:10])}")
        
        # Export comparison as PDF
        st.markdown("---")
        if st.button("📄 Export Comparison Report (PDF)"):
            st.info("PDF export feature coming in Week 4")
        
    except Exception as e:
        st.error(f"Comparison visualization error: {str(e)}")
